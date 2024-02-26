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
    save.image(paste0("workflow/",WILDCARDS[[1]],"_", WILDCARDS[[2]], "_", WILDCARDS[[3]], "_", "create_se_list.RData"))
}
q("no")
print(.libPaths())

# set lib_path to "/usr/local/lib/R/site-library"
.libPaths(c("/usr/local/lib/R/site-library", "/cluster/home/t119797uhn/R/x86_64-pc-linux-gnu-library/4.2"))
library(tximport)
library(GenomicFeatures)
library(rtracklayer)
library(AnnotationDbi)

# load("combine_abundances.RData")
print("Loading Genomic annotation Data")
ANNOTATION_GTF <- INPUT$annotation_file
txdb <- GenomicFeatures::makeTxDbFromGFF(file = ANNOTATION_GTF)
k <- AnnotationDbi::keys(txdb, keytype = "TXNAME")
tx2gene <- AnnotationDbi::select(txdb, k, "GENEID", "TXNAME")
granges <- rtracklayer::import(ANNOTATION_GTF)

print("Loading Kallisto Data")
quant_files <- INPUT$quant_files
quant_files <- quant_files[grepl("h5", quant_files)]

sampleNames <- basename(dirname(quant_files))
names(quant_files) <- sampleNames

quant_files <- quant_files

print("Loading Kallisto Gene Data")
rnaseq.genes <- tximport::tximport(
    quant_files, type = "kallisto", tx2gene = tx2gene, ignoreAfterBar = TRUE, ignoreTxVersion = FALSE)

print("Loading Kallisto Transcript Data")
# countsFromAbundance = c("no", "scaledTPM", "lengthScaledTPM", "dtuScaledTPM"),
countsFromAbundance <- "no"
rnaseq.transcripts <- tximport::tximport(
    quant_files, type = "kallisto", txOut = T, countsFromAbundance = countsFromAbundance)

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
    transcripts_length = rnaseq.transcripts$length)
metadata <- WILDCARDS[which(names(WILDCARDS) != "")]

# 
# sample_metadata <- fread(INPUT$sample_metadata, check.names=T)


# cols <- list(
#     sampleid = "Characteristics.cell.line.",
#     tissueid = "Characteristics.tissue.supergroup.",
#     media = "Characteristics.media.",
#     freeze_media = "Characteristics.freeze.media.",
#     protocol_ref = "Protocol.REF",
#     ENA_ALIAS = "Comment.ENA_ALIAS.",
#     ENA_SAMPLE = "Comment.ENA_SAMPLE.",
#     ENA_EXPERIMENT = "Comment.ENA_EXPERIMENT.",
#     orgnamism_part = "Characteristics.organism.part.",
#     metatastic_site = "Characteristics.metatastic.tissue.",
#     ethnicity = "Characteristics.ethnicity.",
#     assay_name = "Assay.Name",
#     ENA_RUN = "Comment.ENA_RUN.",
#     source = "Performer"
# )

# # rename the columns by replacing the values with the keys
# for (i in seq_along(cols)){
#     data.table::setnames(sample_metadata, cols[[i]], names(cols)[i], skip_absent = TRUE)
# }
# cols <- names(cols)
# sample_metadata <- unique(sample_metadata[, ..cols])



# assay <- rnaseq.genes$abundance
# colnames(assay)[match(sample_metadata$ENA_ALIAS, colnames(assay))] <- sample_metadata$sampleid

rse_list <- lapply(names(assays), function(x){
    assay <- assays[[x]]

    granges <- if(grepl("genes", x)){
        granges[rownames(assay) %in% granges$gene_id & granges$type == "gene",]
    } else {
        granges[rownames(assay) %in% granges$transcript_id & granges$type == "transcript",]
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

print("Saving RDS")
dir.create(dirname(OUTPUT$rse_list), recursive = TRUE, showWarnings = FALSE)
save(rse_list, file = OUTPUT$rse_list, compress = "bzip2", compression_level = 9)


# print("Saving TSVs")
# dir.create(dirname(OUTPUT$assay_tsv[1]), recursive = TRUE, showWarnings = FALSE)
# a_list <- lapply(seq_along(names(assays)), function(x){
#     # write.table(assays[[sub("\\..*", "", basename(x))]], file = x, sep = "\t", quote = FALSE)
#     file_ending <- paste0(names(assays)[x], ".tsv.gz")
#     assay <- assays[[x]]
#     (output_file <- OUTPUT$assay_tsv[grepl(paste0(file_ending, "$"), OUTPUT$assay_tsv)])
#     print(paste0("Writing ", names(assays)[x], " to ", output_file))
#     write.table(
#         assay, 
#         file = output_file,
#         sep = "\t",
#         quote = FALSE
#     )
# })

# write.table(metadata, file = OUTPUT$metadata, sep = "\t", quote = FALSE)

# STORAGE
# $78,330.03
# STANDARD (Cloud Storage)
# $34,535.29

# Total amount of storage
# 1200 TB
# $34,535.29
# Location type
# Region
# N/A
# Location
# Toronto (northamerica-northeast2)
# N/A
# Storage class
# Standard Storage
# N/A
# Source region
# North America
# N/A
# Destination region
# North America
# N/A
# NEARLINE (Cloud Storage)
# $20,771.22

# Data retrieval amount
# 100 TB
# $1,251.28
# Total amount of storage
# 1200 TB
# $19,519.94
# Location type
# Region
# N/A
# Location
# Toronto (northamerica-northeast2)
# N/A
# Storage class
# Nearline Storage
# N/A
# Source region
# North America
# N/A
# Destination region
# North America
# N/A
# COLDLINE (Cloud Storage)
# $13,013.30

# Data retrieval amount
# 100 TB
# $2,502.56
# Total amount of storage
# 1200 TB
# $10,510.74
# Location type
# Region
# N/A
# Location
# Toronto (northamerica-northeast2)
# N/A
# Storage class
# Coldline Storage
# N/A
# Source region
# North America
# N/A
# Destination region
# North America
# N/A
# ARCHIVE (Cloud Storage)
# $10,010.23

# Data retrieval amount
# 100 TB
# $6,256.39
# Total amount of storage
# 1200 TB
# $3,753.84
# Location type
# Region
# N/A
# Location
# Toronto (northamerica-northeast2)
# N/A
# Storage class
# Archive Storage
# N/A
# Source region
# North America
# N/A
# Destination region
# North America
# N/A



# STANDARD:
#     Storage: $34,535.29
#     RETRIEVAL : FREE

# NEARLINE:
#     Storage: $19,519.94
#     RETRIEVAL: $1,251.28

# COLDLINE:
#     Storage: $10,510.74
#     RETRIEVAL: $2,502.56

# ARCHIVE:
#     Storage: $3,753.84
#     RETRIEVAL: $6,256.39
