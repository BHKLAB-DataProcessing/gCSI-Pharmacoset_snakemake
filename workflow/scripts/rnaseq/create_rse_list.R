## ------------------- Parse Snakemake Object ------------------- ##
# This checks if the snakemake object exists. an alternative if running as a script
if(exists("snakemake")){
    INPUT <- snakemake@input
    OUTPUT <- snakemake@output
    WILDCARDS <- snakemake@wildcards
    THREADS <- snakemake@threads

    LOGFILE <- snakemake@log[[1]]
    # this will redirect all output to the logfile
    # split is the equivalent of unix's `tee`
    sink(LOGFILE, type=c("output", "message"), split=TRUE)
    metadata <- WILDCARDS[which(names(WILDCARDS) != "")]
    save.image(paste0("logs/",paste0(metadata, collapse="_"), "_create_se_list.RData"))
    # save.image(paste0("workflow/",WILDCARDS[[1]],"_", WILDCARDS[[2]], "_", WILDCARDS[[3]], "_", "create_se_list.RData"))
}
####################### Generalize functions #######################

#' getTranscripts Function
#'
#' This function imports transcript-level quantification files and performs necessary preprocessing steps.
#'
#' @param quant_files A character vector specifying the paths to the quantification files.
#' @param tx2gene A TxDb object or a character vector specifying the path to the transcript-to-gene mapping file.
#' @param tool A character string specifying the quantification tool used (e.g., "salmon", "kallisto").
#' @param countsFromAbundance A character string specifying whether to derive counts from abundance values ("yes" or "no").
#'
getTranscripts <- function(quant_files, tx2gene, tool, countsFromAbundance = "no"){
    transcripts <- tximport::tximport(
        files = quant_files, 
        type = tool, 
        txOut = TRUE, 
        ignoreTxVersion = FALSE,
        countsFromAbundance=countsFromAbundance
    )
    rownames(transcripts$counts) <- sub("\\|.*", "", rownames(transcripts$counts))
    rownames(transcripts$abundance) <- sub("\\|.*", "", rownames(transcripts$abundance))
    rownames(transcripts$length) <- sub("\\|.*", "", rownames(transcripts$length))

    return(transcripts)
}

#' getGenes Function
#'
#' This function imports gene expression quantification files and performs gene-level summarization using the tximport package.
#'
#' @param quant_files A character vector specifying the paths to the quantification files.
#' @param tx2gene A TxDb object or a named character vector mapping transcript IDs to gene IDs.
#' @param tool A character string specifying the quantification tool used (e.g., "salmon", "kallisto").
#'
getGenes <- function(quant_files, tx2gene, tool){
    print(sprintf("Loading %s Gene Data for %s number of files", tool, length(quant_files)))
    genes <- tximport::tximport(
        quant_files, 
        type = tool, 
        tx2gene = tx2gene, 
        txIn = ifelse(tool == "rsem", FALSE, TRUE), 
        ignoreAfterBar = TRUE, 
        ignoreTxVersion = FALSE
    )
    return(genes)
}


# set lib_path to "/usr/local/lib/R/site-library"
.libPaths(c("/usr/local/lib/R/site-library", "/cluster/home/t119797uhn/R/x86_64-pc-linux-gnu-library/4.2"))
library(tximport)
library(GenomicFeatures)
library(rtracklayer)
library(AnnotationDbi)

print("Loading Genomic annotation Data")
txdb <- GenomicFeatures::makeTxDbFromGFF(file = INPUT$annotation_file)
k <- AnnotationDbi::keys(txdb, keytype = "TXNAME")
tx2gene <- AnnotationDbi::select(txdb, k, "GENEID", "TXNAME")
GRanges <- rtracklayer::import(INPUT$annotation_file)


# Tool should be one of "salmon", "kallisto"
tool <- WILDCARDS$TOOL
stopifnot(tool %in% c("salmon", "kallisto", "rsem"))

quant_files <- INPUT$quant_files
names(quant_files) <- basename(dirname(quant_files))

if(tool %in% c("salmon", "kallisto")){
    rnaseq.genes<- getGenes(quant_files, tx2gene, tool)
    rnaseq.transcripts <- getTranscripts(quant_files, tx2gene, tool, countsFromAbundance = "no")

}else if(tool == "rsem"){
    # split into genes and transcripts
    rsem_genes <- quant_files[grepl("genes.results", quant_files)]
    rsem_transcripts <- quant_files[grepl("isoforms.results", quant_files)]

    rnaseq.genes <- getGenes(rsem_genes, tx2gene, tool)
    rnaseq.transcripts <- getTranscripts(rsem_transcripts, tx2gene, tool, countsFromAbundance = "no")
}


### Save the data
assays <- list(
    genes = rnaseq.genes$abundance,
    genes_counts = rnaseq.genes$counts,
    genes_length = rnaseq.genes$length,
    transcripts = rnaseq.transcripts$abundance,
    transcripts_counts = rnaseq.transcripts$counts,
    transcripts_length = rnaseq.transcripts$length
)
str(assays)
metadata <- WILDCARDS[which(names(WILDCARDS) != "")]



rse_list <- lapply(names(assays), function(x){
    assay <- assays[[x]]

    if(grepl("genes", x)){
        granges <-  GRanges[GRanges$gene_id %in% rownames(assay) & GRanges$type == "gene",]
        assay <- assay[rownames(assay) %in% granges$gene_id,]
    } else {
        granges <-  GRanges[GRanges$transcript_id %in% rownames(assay) & GRanges$type == "transcript",]
        assay <- assay[rownames(assay) %in% granges$transcript_id,]
    }

    SummarizedExperiment::SummarizedExperiment(
        assays = list(
            expr = assay
        ),
        rowRanges = granges,
        metadata = metadata
    )
})
names(rse_list) <- names(assays)
print(rse_list)


print(paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " Saving RData to ", OUTPUT$rse_list))
object.size(rse_list)  |> print()

dir.create(dirname(OUTPUT$rse_list), recursive = TRUE, showWarnings = FALSE)
save(rse_list, file = OUTPUT$rse_list, compress = "bzip2", compression_level = 9)


# quant_files <- INPUT$quant_files[1:10]
# names(quant_files) <- basename(dirname(quant_files))
# rnaseq.genes<- getGenes(quant_files, tx2gene, tool)

# rnaseq.transcripts <- getTranscripts(quant_files, tx2gene, tool, countsFromAbundance = "no")

# # print(sprintf("Loading %s Gene Data", tool))
# rnaseq.genes <- tximport::tximport(
#     quant_files, 
#     type = tool,
#     tx2gene = tx2gene, 
#     ignoreAfterBar = TRUE, 
#     ignoreTxVersion = FALSE)
    

# print(sprintf("Loading %s Transcript Data", tool))
# countsFromAbundance <- "no"
# rnaseq.transcripts <- tximport::tximport(
#     quant_files, 
#     type = tool, 
#     txOut = TRUE,
#     countsFromAbundance = countsFromAbundance)

# print("Configuring transcripts Rownames")
# rownames(rnaseq.transcripts$counts) <- sub("\\|.*", "", rownames(rnaseq.transcripts$counts))
# rownames(rnaseq.transcripts$abundance) <- sub("\\|.*", "", rownames(rnaseq.transcripts$abundance))
# rownames(rnaseq.transcripts$length) <- sub("\\|.*", "", rownames(rnaseq.transcripts$length))