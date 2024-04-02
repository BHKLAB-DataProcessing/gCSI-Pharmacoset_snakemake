## ------------------- Parse Snakemake Object ------------------- ##
# Check if the "snakemake" object exists
if(exists("snakemake")){
    INPUT <- snakemake@input
    OUTPUT <- snakemake@output
    WILDCARDS <- snakemake@wildcards
    THREADS <- snakemake@threads

    # setup logger if log file is provided
    if(length(snakemake@log)>0) 
        sink(
            file = snakemake@log[[1]], 
            append = FALSE, 
            type = c("output", "message"), 
            split = TRUE
    )

    file.path("resources", paste0(snakemake@rule, ".RData")) |> 
        save.image()
}else{
    file.path("resources", "annotate_treatmentMetadata.RData") |>
        load()
}

#############################################################################
# Load INPUT
#############################################################################
rawTreatmentMetadata <- unique(data.table::fread(INPUT$rawTreatmentMetadata, header = TRUE, sep = "\t"))
rawTreatmentMetadata


#############################################################################
# Main Script
(compounds_to_cids <- 
  rawTreatmentMetadata[, 
    AnnotationGx::mapCompound2CID(
        names = Norm_DrugName,
        first = TRUE
        )
      ]
)


properties <- c("Title", "MolecularFormula", "InChIKey", "MolecularWeight")

cids_to_properties <- compounds_to_cids[
  !is.na(cids), 
  {
    AnnotationGx::mapCID2Properties(
      ids = cids,
      properties = properties
      )
  }
]


gCSI_treatmentMetadata_annotated <- merge(
  cids_to_properties,
  compounds_to_cids,
  by.y = "cids",
  by.x = "CID",
  all = TRUE
)
gCSI_treatmentMetadata_annotated


# Renaming Columns
# ----------------

old_names <- c("DrugName", "Norm_DrugName")
new_names <- c("gCSI.treatmentid", "gCSI.norm_treatmentid")
data.table::setnames(rawTreatmentMetadata, old = old_names, new = new_names)

names(gCSI_treatmentMetadata_annotated) <- paste("pubchem", names(gCSI_treatmentMetadata_annotated), sep = ".")

gCSI_treatmentMetadata <- merge(
    rawTreatmentMetadata,
    gCSI_treatmentMetadata_annotated,
    by.x = "gCSI.norm_treatmentid",
    by.y = "pubchem.name",
    all.x = TRUE
)

# Annotate with Unichem
# ----------------------
str(gCSI_treatmentMetadata)

unichem_sources <- AnnotationGx::getUnichemSources(T)
data.table::setkey(unichem_sources, Name)
sources_of_interest <- c("chembl", "drugbank", "chebi", "phamgkb", "lincs", "clinicaltrials", "nih_ncc", "fdasrs", "pharmgkb", "rxnorm")
sourceID <- unichem_sources[Name == "pubchem", SourceID]



message("\n\nAnnotating with unichem...")
annotations <- lapply(gCSI_treatmentMetadata$pubchem.CID, function(x){
  tryCatch({
    result <- AnnotationGx::queryUnichemCompound(type = "sourceID", compound = x, sourceID = sourceID)

    subset <- result$External_Mappings[Name %in% sources_of_interest, .(compoundID, Name)]
    # make Name the column names and the values the compoundID 
    subset$cid <- x
    dcast(subset, cid ~ Name, value.var = "compoundID", fun.aggregate = list)
  }, error = function(e) NULL)
  } 
  ) |> data.table::rbindlist(fill = T)
show(annotations)



unichem_mappings <- copy(annotations)
# for each column, if its a list then make it a string with a comma separator
for(col in names(unichem_mappings)){
  if(is.list(unichem_mappings[[col]])){
    unichem_mappings[[col]] <- sapply(unichem_mappings[[col]], function(x) paste(x, collapse = ","))
  }
}
# Rename columns like drugbank to unichem.DrugBank etc, using the unichem_sources "NameLabel" column 
names(unichem_mappings) <- paste("unichem", unichem_sources[names(unichem_mappings), gsub(" ", "_", NameLabel)], sep = ".")


all_annotated_treatmentMetadata <- merge(
    gCSI_treatmentMetadata, 
    unichem_mappings, 
    by.x = "pubchem.CID", 
    by.y = "unichem.NA",
    all.x = T
)

str(all_annotated_treatmentMetadata)

# Annotate with ChEMBL
# ---------------------
message("\n\nAnnotating with ChEMBL using Unichem-obtained ChEMBL IDs")
chembl_mechanisms_dt <- all_annotated_treatmentMetadata[, AnnotationGx::getChemblMechanism(unichem.ChEMBL)]

chembl_cols_of_interest <- c(
        "molecule_chembl_id",  "parent_molecule_chembl_id", "target_chembl_id", "record_id", 
        "mechanism_of_action", "mechanism_comment", "action_type"
    )

all_annotated_treatmentMetadata <- merge(
    all_annotated_treatmentMetadata, 
    chembl_mechanisms_dt[!duplicated(molecule_chembl_id), ..chembl_cols_of_interest], 
    by.x = "unichem.ChEMBL",
    by.y = "molecule_chembl_id", 
    all.x = TRUE
)

data.table::setnames(
    all_annotated_treatmentMetadata, 
    chembl_cols_of_interest, 
    paste0("chembl.", chembl_cols_of_interest), 
    skip_absent = TRUE
)

all_annotated_treatmentMetadata <- all_annotated_treatmentMetadata[!duplicated(pubchem.CID),]

str(all_annotated_treatmentMetadata)

## --------------------- Save OUTPUT ------------------- ##
outputfile <- OUTPUT$treatmentMetadata

message("Saving treatmentMetadata to: ", outputfile)
dir.create(dirname(outputfile), showWarnings = FALSE, recursive = TRUE)

data.table::fwrite(
    x = all_annotated_treatmentMetadata,
    file = outputfile,
    quote = TRUE,
    sep = "\t",
    na = "NA",
    col.names = TRUE
)
