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
    file.path("resources", "create_MultiAssayExperiment.RData") |>
        load()
}

library(MultiAssayExperiment)
#############################################################################
# Load INPUT
#############################################################################

# check if the rse_list variable exists 
if(!exists("rse_list")){
     load(INPUT$rnaseq) 
}
show(rse_list)

sampleMetadata <- data.table::fread(INPUT$sampleMetadata, header = TRUE, sep = "\t", check.names = T)
str(sampleMetadata)

#############################################################################
# Main Script

(mae_colData <- sampleMetadata[, .(
    gCSI.sampleid,
    EGA.ENA_ALIAS,
    cellosaurus.cellLineName,
    cellosaurus.accession
)])

allColnamesInColData <- sapply(rse_list, function(se){
    all(colnames(se) %in% mae_colData$EGA.ENA_ALIAS)
})

if(!all(allColnamesInColData)){
    stop("Not all the sample names in the RSE list are present in the colData")
}


renamed_rse_list <- lapply(rse_list, function(se){
    colnames(se) <- mae_colData$gCSI.sampleid[match(colnames(se), mae_colData$EGA.ENA_ALIAS)]

    assay_temp <- assay(se, "expr")
    # order the columns in the assay matrix alphabetically

    assay_temp <- assay_temp[, order(colnames(assay_temp))]

    # reassign the assay matrix
    new_colData <- mae_colData[match(colnames(assay_temp), mae_colData$gCSI.sampleid),]
    new_colData[, sampleid := gCSI.sampleid ]
    new_colData[, batchid := NA_character_]
    
    assay(se, "expr", withDimnames= FALSE) <- assay_temp
    new_se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(
            expr = assay_temp
        ),
        rowRanges = rowRanges(se),
        colData = new_colData,
        metadata = metadata(se)
    )
    return(new_se)
}) 
rnaseq_se <- c("genes", "transcripts", "genes_counts", "transcripts_counts", "genes_length", "transcripts_length")

renamed_rse_list <- sapply(rnaseq_se, function(se_name){
    message(paste("Checking", se_name))
    se <- renamed_rse_list[[se_name]]
    metadata(se)$annotation <- "rnaseq"
    se
})

renamed_rse_list

# Build MultiAssayExperiment
# --------------------------

ExpList <- MultiAssayExperiment::ExperimentList(renamed_rse_list) 
message(paste("ExperimentList:", capture.output(show(ExpList)), sep = "\n\t"))

mae_colData[, sampleid := gCSI.sampleid]
mae_colData[, batchid := 1]

data.table::setkeyv(mae_colData, "sampleid")

sampleMapList <- lapply(renamed_rse_list, function(se){
    data.frame(
        primary = colnames(se),
        colname = colnames(se),
        stringsAsFactors = FALSE
    )
})

names(sampleMapList) <- names(ExpList)
message(paste("Sample map list:", capture.output(str(sampleMapList)), sep = "\n\t"))

# Metadata List
# go through each experiment, extract the metadata and add it to a list
metadata_list <- lapply(renamed_rse_list, function(se){
    metadata_ <- slot(se, "metadata")
    metadata_
})
metadata_list



colData <- as.data.frame(
    mae_colData,
    row.names = mae_colData$gCSI.sampleid
)

message(sprintf("Column data has %d rows and %d columns", nrow(colData), ncol(colData)))
str(colData)

mae <- MultiAssayExperiment::MultiAssayExperiment(
    experiments = ExpList,
    colData = colData,
    sampleMap = MultiAssayExperiment::listToMap(sampleMapList),
    metadata = metadata_list
)

## --------------------- Save OUTPUT ------------------- ##

message("Saving MultiAssayExperiment to: ", OUTPUT$mae)
dir.create(
    path = dirname(OUTPUT$mae),
    showWarnings = FALSE,
    recursive = TRUE
)

saveRDS(mae, file = OUTPUT$mae)
