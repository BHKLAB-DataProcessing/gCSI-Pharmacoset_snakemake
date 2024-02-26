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
    save.image(paste0("workflow/",WILDCARDS[[1]],"_", WILDCARDS[[2]], "_","combine_salmon_quant.RData"))
}

print(.libPaths())

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

print("Loading Salmon Gene Data")
quant_files <- INPUT$salmon_quant
names(quant_files) <- basename(dirname(quant_files))
rnaseq.genes <- tximport::tximport(
    quant_files, 
    type = "salmon", 
    tx2gene = tx2gene, 
    ignoreAfterBar = TRUE, 
    ignoreTxVersion = FALSE)

print("Loading Salmon Transcript Data")
countsFromAbundance <- "no"
rnaseq.transcripts <- tximport::tximport(
    quant_files, 
    type = "salmon", 
    txOut = TRUE,
    countsFromAbundance = countsFromAbundance)

print("Configuring transcripts Rownames")
rownames(rnaseq.transcripts$counts) <- sub("\\|.*", "", rownames(rnaseq.transcripts$counts))
rownames(rnaseq.transcripts$abundance) <- sub("\\|.*", "", rownames(rnaseq.transcripts$abundance))
rownames(rnaseq.transcripts$length) <- sub("\\|.*", "", rownames(rnaseq.transcripts$length))

### Save the data
assays <- list(
    genes = rnaseq.genes$abundance,
    genes_counts = rnaseq.genes$counts,
    genes_length = rnaseq.genes$length,
    transcripts = rnaseq.transcripts$abundance,
    transcripts_counts = rnaseq.transcripts$counts,
    transcripts_length = rnaseq.transcripts$length
)

metadata <- WILDCARDS[which(names(WILDCARDS) != "")]

output_RDS <- list(
    tximport = list(
        rnaseq.genes = rnaseq.genes,
        rnaseq.transcripts = rnaseq.transcripts
    ),
    GRanges = granges,
    assays = assays,
    metadata = metadata
)

rnaseq.genes <- output_RDS$tximport$rnaseq.genes
rnaseq.transcripts <- output_RDS$tximport$rnaseq.transcripts
GRanges <- output_RDS$GRanges

rse <- SummarizedExperiment::SummarizedExperiment(
    assays = list(
        genes = rnaseq.genes$abundance
    ),
    rowRanges = GRanges_genes
)



print("Saving RDS")
dir.create(dirname(OUTPUT$RDS), recursive = TRUE, showWarnings = FALSE)
saveRDS(output_RDS, file = OUTPUT$RDS)


print("Saving TSVs")
dir.create(dirname(OUTPUT$assay_tsv[1]), recursive = TRUE, showWarnings = FALSE)
a_list <- lapply(seq_along(names(assays)), function(x){
    # write.table(assays[[sub("\\..*", "", basename(x))]], file = x, sep = "\t", quote = FALSE)
    file_ending <- paste0(names(assays)[x], ".tsv.gz")
    assay <- assays[[x]]
    (output_file <- OUTPUT$assay_tsv[grepl(paste0(file_ending, "$"), OUTPUT$assay_tsv)])
    print(paste0("Writing ", names(assays)[x], " to ", output_file))
    write.table(
        assay, 
        file = output_file,
        sep = "\t",
        quote = FALSE,
        header = TRUE,
    )
})


write.table(metadata, file = OUTPUT$metadata, sep = "\t", quote = FALSE, col.names = TRUE)

#### RSEM


rnaseq.genes <- tximport::tximport(
    genes, 
    type = "rsem", 
    tx2gene = tx2gene, 
    ignoreAfterBar = TRUE, 
    ignoreTxVersion = FALSE,
    txIn = FALSE, 
    txOut = FALSE
)


rnaseq.transcripts <- tximport::tximport(
    transcripts, 
    type = "rsem", 
    tx2gene = tx2gene, 
    ignoreAfterBar = TRUE, 
    ignoreTxVersion = FALSE,
    txIn = TRUE, 
    txOut = TRUE
)

print("Configuring transcripts Rownames")
rownames(rnaseq.transcripts$counts) <- sub("\\|.*", "", rownames(rnaseq.transcripts$counts))
rownames(rnaseq.transcripts$abundance) <- sub("\\|.*", "", rownames(rnaseq.transcripts$abundance))
rownames(rnaseq.transcripts$length) <- sub("\\|.*", "", rownames(rnaseq.transcripts$length))


### Save the data
assays <- list(
    genes = rnaseq.genes$abundance,
    genes_counts = rnaseq.genes$counts,
    genes_length = rnaseq.genes$length,
    transcripts = rnaseq.transcripts$abundance,
    transcripts_counts = rnaseq.transcripts$counts,
    transcripts_length = rnaseq.transcripts$length
)
