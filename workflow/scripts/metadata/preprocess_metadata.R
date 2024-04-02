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
    file.path("resources", "preprocess_metadata.RData") |>
        load()
}


#############################################################################
# Load INPUT
#############################################################################
tarGZ_file <- INPUT$rawData

message("Untar the file and list files")
# untar the file and list files
untar(tarGZ_file, exdir = "rawdata/tmp") 

# list files
files <- list.files("rawdata/tmp", full.names = TRUE, recursive = TRUE)
metricsFile <- grep("GRmetrics", files, value = TRUE) 
rawData <- data.table::fread(metricsFile, header = TRUE, sep = "\t")

show(rawData)

#############################################################################
# Main Script
# Subset the metadata

sampleMetadata <- rawData[, c(
    "CellLineName", "PrimaryTissue", "Norm_CellLineName"
)]

treatmentMetadata <- rawData[, c(
    "DrugName", "Norm_DrugName"
)]

#############################################################################
# SAVE OUTPUT
#############################################################################

rawSampleMetadata <- OUTPUT$rawSampleMetadata
message("Saving raw sample metadata to ", rawSampleMetadata)
data.table::fwrite(sampleMetadata, rawSampleMetadata, sep = "\t", quote = FALSE)

rawTreatmentMetadata <- OUTPUT$rawTreatmentMetadata
message("Saving raw treatment metadata to ", rawTreatmentMetadata)
data.table::fwrite(treatmentMetadata, rawTreatmentMetadata, sep = "\t", quote = FALSE)