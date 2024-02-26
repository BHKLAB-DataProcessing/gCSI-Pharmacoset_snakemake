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
}else {
    # if not running as a snakemake script, then set the variables manually
# assay_tsv = multiext(
#             "results/data/rnaseq/salmon_v{SALMON_version}_{ref_build}.{ref_version}/assay",
#             "-genes.tsv.gz",
#             "-genes_counts.tsv.gz",
#             "-genes_length.tsv.gz",
#             "-transcripts.tsv.gz",
#             "-transcripts_counts.tsv.gz",
#             "-transcripts_length.tsv.gz",
#         ),
    SALMON_version <- "1.8.0"
    KALLISTO_version <- "0.46.1"
    ref_build <- "GRCh38"
    ref_version <- "45"
    salmon_prefix <- sprintf("results/data/rnaseq/salmon_v%s_%s.%s/", SALMON_version, ref_build, ref_version)
    kallisto_prefix <- sprintf("results/data/rnaseq/kallisto_v%s_%s.%s/", KALLISTO_version, ref_build, ref_version)
    assay_tsvs <- paste0(kallisto_prefix, c(
        "assay-genes.tsv.gz",
        "assay-genes_counts.tsv.gz",
        "assay-genes_length.tsv.gz",
        "assay-transcripts.tsv.gz",
        "assay-transcripts_counts.tsv.gz",
        "assay-transcripts_length.tsv.gz"
    ))
}

dt <- data.table::fread(assay_tsvs[1], header = TRUE, sep = "\t")

    # assays <- INPUT$assay_tsv