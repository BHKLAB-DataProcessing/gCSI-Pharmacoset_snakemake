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
    file.path("resources", "build_PharmacoSet.RData") |>
        load()
}



#############################################################################
# Load INPUT
#############################################################################
message("Loading: ", INPUT$treatmentMetadata)
treatmentMetadata <- data.table::fread(INPUT$treatmentMetadata, check.names = T, sep = "\t", header = T)

message("Loading: ", INPUT$sampleMetadata)
sampleMetadata <- data.table::fread(INPUT$sampleMetadata, check.names = T, sep = "\t", header = T)

message("Loading: ", INPUT$mae)
mae <- readRDS(INPUT$mae)

message("Loading: ", INPUT$tre)
tre <- readRDS(INPUT$tre)


#############################################################################
# Main Script
data.table::setkeyv(sampleMetadata, "gCSI.sampleid")
sampleMetadata[, sampleid := gCSI.sampleid]
sampleMetadata <- sampleMetadata[!duplicated(sampleid),][!is.na(sampleid),]

sample <- as.data.frame(
    sampleMetadata, 
    row.names = sampleMetadata[, sampleid]
)


data.table::setkeyv(treatmentMetadata, "gCSI.treatmentid")
treatmentMetadata[, treatmentid := gCSI.treatmentid]

treatment <- as.data.frame(
    treatmentMetadata, 
    row.names = treatmentMetadata[, treatmentid]
)


name <- "gCSI"
pset <- PharmacoGx::PharmacoSet2(
    name = name,
    treatment = treatment,
    sample = sample,
    molecularProfiles = mae,
    treatmentResponse = tre,
    perturbation = list(),
    curation = list(
        sample = sample, 
        treatment = treatment, 
        tissue = data.frame()),
    datasetType = "sensitivity"
)
message(paste(capture.output(show(pset)), collapse = "\n\t"))


## ----------------------------------------------------- ##




## --------------------- Save OUTPUT ------------------- ##

message("Object Size (pset):")
object.size(pset) |> print(unit = "auto")
print(paste("Saving PharmacoSet object to", OUTPUT[[1]]))
dir.create(dirname(OUTPUT[[1]]), recursive = TRUE, showWarnings = FALSE)
saveRDS(pset, file = OUTPUT[[1]])