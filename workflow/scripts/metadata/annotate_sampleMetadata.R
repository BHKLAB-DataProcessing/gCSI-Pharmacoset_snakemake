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
    file.path("resources", "annotate_sampleMetadata.RData") |>
        load()
}

options("mc.cores" = THREADS)
options("log_level" = "INFO")   # AnnotationGx logging level

#############################################################################
# LOAD INPUT
#############################################################################
rawSampleMetadata <- unique(data.table::fread(INPUT$rawSampleMetadata, header = TRUE, sep = "\t"))
rawSampleMetadata <- rawSampleMetadata[order(CellLineName)]
head(rawSampleMetadata)

EGA_sampleMetadata <- data.table::fread(INPUT$EGA_sampleMetadata, header = TRUE, sep = "\t", check.names = T)
head(EGA_sampleMetadata)

#############################################################################
# MAIN SCRIPT
#############################################################################

# PREPROCESSSING ALL THE SAMPLE METADATA
# -------------------------------------

# Add the "gCSI" prefix to the column names
names(rawSampleMetadata) <- paste("gCSI", names(rawSampleMetadata), sep = ".")
message("\n\nNumber of rows in the raw sample metadata: ", nrow(rawSampleMetadata))

# Add the "EGA" prefix to the column names
## ----------------------------------------------------- ##
EGA_sampleMetadata_subset <- EGA_sampleMetadata[, .(
    Characteristics.cell.line., Comment.ENA_ALIAS., 
    Characteristics.tissue.supergroup., Characteristics.metatastic.tissue., Factor.Value.disease.,
    Characteristics.media., Characteristics.freeze.media., 
    Characteristics.age., Characteristics.sex., Characteristics.ethnicity., Protocol.REF,
    Assay.Name, Scan.Name
)]
EGA_sampleMetadata_subset
old_names <- c(
    "Characteristics.cell.line.", "Comment.ENA_ALIAS.",
    "Characteristics.tissue.supergroup.", "Characteristics.metatastic.tissue.", "Factor.Value.disease.",
    "Characteristics.media.", "Characteristics.freeze.media.",
    "Characteristics.age.", "Characteristics.sex.", "Characteristics.ethnicity.", "Protocol.REF",
    "Assay.Name", "Scan.Name"
)

new_names <- c(
    "CellLineName", "ENA_ALIAS",
    "TissueSupergroup", "MetastaticTissue", "Disease",
    "Media", "FreezeMedia",
    "Age", "Sex", "Ethnicity", "Protocol_REF",
    "AssayName", "ScanName"
)

data.table::setnames(EGA_sampleMetadata_subset, old = old_names, new = new_names)

names(EGA_sampleMetadata_subset) <- paste("EGA", names(EGA_sampleMetadata_subset), sep = ".")

message("\n\nNumber of rows in the EGA sample metadata: ", nrow(EGA_sampleMetadata))
show(EGA_sampleMetadata_subset)


# MERGE THE SAMPLE METADATA
# -------------------------
message("\n\nMerging the sample metadata")
gCSI_sampleMetadata <- merge(
    x = rawSampleMetadata, 
    y = EGA_sampleMetadata_subset, 
    by.x = "gCSI.CellLineName", 
    by.y = "EGA.CellLineName", 
    all = TRUE
)[!duplicated(gCSI.CellLineName),]
show(gCSI_sampleMetadata)


# Annotate with AnnotationGx
# --------------------------
message("\n\n Getting Accessions with cellosaurus...")
sample_accessions <- AnnotationGx::mapCell2Accession(
    ids = gCSI_sampleMetadata$gCSI.CellLineName 
)
data.table::setnames(sample_accessions, "query", "gCSI.sampleid")
str(sample_accessions)
print("\n\nMissing samples in gCSI:")

name <- sample_accessions[is.na(accession),unique(gCSI.sampleid)]
result <- AnnotationGx::mapCell2Accession(name, fuzzy = T)
data.table::setnames(result, "query", "gCSI.sampleid")
sample_accessions <- data.table::rbindlist(
  list(
    unique(sample_accessions[!is.na(accession)]),
    result
  )
) |> unique()
sample_accessions <- sample_accessions[order(gCSI.sampleid)]

message("\n\nNumber of rows in the sample accessions: ", nrow(sample_accessions))

message("\n\nAnnotating with cellosaurus...")
annotated_accessions <- sample_accessions[!is.na(accession), AnnotationGx::annotateCellAccession(accession)]
annotated_accessions[, synonyms := sapply(synonyms, function(x) paste(x, collapse = "; "))]
annotated_accessions[, diseases := sapply(diseases, function(x) paste(x, collapse = "; "))]
annotated_accessions[, c("crossReferences", "hierarchy", "comments") := NULL]
names(annotated_accessions) <- paste0("cellosaurus.", names(annotated_accessions))
annotated_accessions <- unique(annotated_accessions)

cellosaurus_annotations <- merge(
    annotated_accessions,
    sample_accessions,
    by.x = c("cellosaurus.cellLineName", "cellosaurus.accession"),
    by.y = c("cellLineName", "accession"),
)


final_annotated <- 
    merge(
        x = gCSI_sampleMetadata, 
        y = cellosaurus_annotations,
        by.x = c("gCSI.CellLineName"), 
        by.y = c("gCSI.sampleid"),
        all.x = TRUE
    ) |> unique()

final_annotated[, gCSI.sampleid := gCSI.CellLineName]
final_annotated[, sampleid := gCSI.sampleid]
final_annotated


#############################################################################
# SAVE OUTPUT
#############################################################################
outputfile <- OUTPUT$sampleMetadata

message("Saving sampleMetadata to: ", outputfile)

dir.create(
    path = dirname(outputfile),
    showWarnings = FALSE,
    recursive = TRUE
)

data.table::fwrite(
    x = final_annotated,
    file = outputfile,
    quote = TRUE,
    sep = "\t",
    na = "NA",
    row.names = FALSE
)

