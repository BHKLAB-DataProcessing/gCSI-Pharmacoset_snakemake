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
    file.path("resources", "create_treatmentResponseExperiment.RData") |>
        load()
}


#############################################################################
# Load INPUT
#############################################################################

treatmentMetadata <- data.table::fread(INPUT$treatmentMetadata, check.names = T, sep = "\t", header = T)

sampleMetadata <- data.table::fread(INPUT$sampleMetadata, check.names = T, sep = "\t", header = T)


tsv_path <- INPUT$tsv_path
tr_dir <- dirname(tsv_path)
unzip_dir <- file.path(tr_dir, "extracted")

if(!dir.exists(unzip_dir)) dir.create(unzip_dir)

message("Extracting TSV")
system(paste("tar -xzf", tsv_path, "--strip-components 1 -C", unzip_dir))
tsv_files <- list.files(unzip_dir, pattern = "*.tsv", full.names = TRUE)
input_dt_list <- lapply(tsv_files, data.table::fread)
names(input_dt_list) <- basename(tsv_files)

GRValues <- input_dt_list[["gCSI_GRvalues_v1.3.tsv"]]
#############################################################################
# Main Script



raw_dt <- GRValues[,
    .(
    gCSI.treatmentid = Norm_DrugName,
    gCSI.sampleid = Norm_CellLineName,
    ExperimentNumber,
    dose = ((10^log10Concentration)/1e-6),
    viability = relative_cell_count * 100,
    # expid = paste(Norm_CellLineName, Norm_DrugName, ExperimentNumber, sep = "_"),
    GR,
    activities.mean,
    activities.std
)]
treatmentNames <- raw_dt$gCSI.treatmentid |> unique() |> sort()
message("\n\n Treatment Names:\n\t", paste(treatmentNames, collapse = "\n\t"))

message("\n\n Treatment Metadata:")
str(treatmentMetadata)

cols_ <- c("gCSI.norm_treatmentid", "gCSI.treatmentid")
tre_treatmentMetadata <- treatmentMetadata[, ..cols_]

missing <- setdiff(treatmentNames, tre_treatmentMetadata$gCSI.norm_treatmentid)

if(length(missing) > 0){
    message("\n\n Missing Treatments:")
    message("\t", paste(missing, collapse = "\n\t"))}

message("\n\n Removing missing treatments from the raw data")
raw_dt <- raw_dt[!gCSI.treatmentid %in% missing]

raw_dt[, 
    treatmentid := tre_treatmentMetadata[match(raw_dt$gCSI.treatmentid, tre_treatmentMetadata$gCSI.norm_treatmentid), gCSI.treatmentid]
]

## SAMPLES

sampleNames <- raw_dt$gCSI.sampleid |> unique() |> sort()

missingSamples <- setdiff(sampleNames, sampleMetadata$gCSI.Norm_CellLineName) 

if(length(missingSamples) > 0){
    message("\n\n Missing Samples:")
    message("\t", paste(missingSamples, collapse = "\n\t"))
}

raw_dt[
    , sampleid :=  sampleMetadata[match(raw_dt$gCSI.sampleid, sampleMetadata$gCSI.Norm_CellLineName), gCSI.sampleid]
]

#### DOUBLE CHECKING
# assert no missing values in treatmentid and sampleid
stopifnot(all(!is.na(raw_dt$treatmentid)))
stopifnot(all(!is.na(raw_dt$sampleid)))


message(paste0("Loading TREDataMapper"))
TREDataMapper <- CoreGx::TREDataMapper(rawdata=raw_dt)

CoreGx::rowDataMap(TREDataMapper) <- list(
    id_columns = (c("treatmentid", "gCSI.treatmentid", "ExperimentNumber", "dose")),
    mapped_columns = c())

CoreGx::colDataMap(TREDataMapper) <- list(
    id_columns = c("sampleid", "gCSI.sampleid"),
    mapped_columns = c())

CoreGx::assayMap(TREDataMapper) <- list(
    sensitivity = list(
        c("treatmentid", "sampleid", "gCSI.treatmentid", "gCSI.sampleid", "ExperimentNumber", "dose"),
        c("viability")),
    summary = list(
        c("treatmentid", "sampleid", "gCSI.treatmentid", "gCSI.sampleid", "ExperimentNumber", "dose"),
        c("GR", "activities.mean", "activities.std")))

gCSI_tre <- CoreGx::metaConstruct(TREDataMapper)
gCSI_tre

message("Computing Curve Fits...")
tre_fit <- gCSI_tre |> CoreGx::endoaggregate(
    {  # the entire code block is evaluated for each group in our group by
        # 1. fit a log logistic curve over the dose range
        fit <- PharmacoGx::logLogisticRegression(dose, viability,
            viability_as_pct=TRUE)
        # 2. compute curve summary metrics
        ic50 <- PharmacoGx::computeIC50(dose, Hill_fit=fit)
        aac <- PharmacoGx::computeAUC(dose, Hill_fit=fit)
        # 3. assemble the results into a list, each item will become a
        #   column in the target assay.
        list(
            HS=fit[["HS"]],
            E_inf = fit[["E_inf"]],
            EC50 = fit[["EC50"]],
            Rsq=as.numeric(unlist(attributes(fit))),
            aac_recomputed=aac,
            ic50_recomputed=ic50
        )
    },
    assay="sensitivity",
    target="profiles",
    enlist=FALSE,  # this option enables the use of a code block for aggregation
    by=c("treatmentid", "sampleid"),
    nthread=THREADS  # parallelize over multiple cores to speed up the computation
)
show(tre_fit)

# Reading in Metrics:
# message("Reading in Metrics...")

# GRMetrics <- input_dt_list[["gCSI_GRmetrics_v1.3.tsv"]]
# # names(GRMetrics)
# # names(GRMetrics)
# #  [1] "DrugName"              "CellLineName"          "PrimaryTissue"
# #  [4] "ExperimentNumber"      "TrtDuration"           "doublingtime"
# #  [7] "maxlog10Concentration" "GR_AOC"                "meanviability"
# # [10] "GRmax"                 "Emax"                  "GRinf"
# # [13] "GEC50"                 "GR50"                  "R_square_GR"
# # [16] "pval_GR"               "flat_fit_GR"           "h_GR"
# # [19] "N_conc"                "log10_conc_step"       "Norm_CellLineName"
# # [22] "Norm_DrugName"         "GR_05uM_fit"

# published_profiles <- unique(
#     GRMetrics[,
#     .(
#     gCSI.treatmentid = Norm_DrugName,
#     gCSI.sampleid = Norm_CellLineName,
#     mean_viability = meanviability,
#     GR_AOC,
#     GRmax,
#     Emax,
#     GRinf,
#     GEC50,
#     GR50)]
# )
# published_profiles[
#     , sampleid := sampleMetadata[match(published_profiles$gCSI.sampleid, sampleMetadata$gCSI.Norm_CellLineName), gCSI.sampleid]
# ]
# published_profiles[
#     , treatmentid := tre_treatmentMetadata[match(published_profiles$gCSI.treatmentid, tre_treatmentMetadata$gCSI.norm_treatmentid), gCSI.treatmentid]
# ]

# CoreGx::assay(tre_fit, "profiles_published") <- published_profiles
# published_profiles <- unique(
#     GRMetrics[,
#     .(
#     gCSI.treatmentid = Norm_DrugName,
#     gCSI.sampleid = Norm_CellLineName,
#     mean_viability = meanviability,
#     GR_AOC,
#     GRmax,
#     Emax,
#     GRinf,
#     GEC50,
#     GR50)])

# CoreGx::assay(gCSI_tre, "profiles_published") <- published_profiles




## ----------------------------------------------------- ##




## --------------------- Save OUTPUT ------------------- ##

outputfile <- OUTPUT$tre
dir.create(
    path = dirname(outputfile),
    showWarnings = FALSE,
    recursive = TRUE
)

message("Saving TRE to: ", outputfile)

saveRDS(gCSI_tre, outputfile)
