
## ------------------- Parse Snakemake Object ------------------- ##
if(exists("snakemake")){
    INPUT <- snakemake@input
    OUTPUT <- snakemake@output
    WILDCARDS <- snakemake@wildcards
    THREADS <- snakemake@threads
    tsv_path <- INPUT$metadata_tsv
    rds_path <- INPUT$metadata_rds
}else{
    tsv_path <- "rawdata/TreatmentResponse_and_Metadata/gCSI_GRdata_v1.3.tsv.tar.gz"
    rds_path <- "rawdata/TreatmentResponse_and_Metadata/gCSI_GRdata_v1.3.rds.tar.gz"
}
save.image()
library(data.table)

tr_dir <- dirname(tsv_path)
unzip_dir <- file.path(tr_dir, "extracted")

# ------------------- Extract TSV ------------------- #
if(!dir.exists(unzip_dir)) dir.create(unzip_dir)

message("Extracting TSV")
system(paste("tar -xzf", tsv_path, "--strip-components 1 -C", unzip_dir))
tsv_files <- list.files(unzip_dir, pattern = "*.tsv", full.names = TRUE)
input_dt_list <- lapply(tsv_files, data.table::fread)
names(input_dt_list) <- basename(tsv_files)


GRValues <- input_dt_list[["gCSI_GRvalues_v1.3.tsv"]]

# names(GRValues)
#  [1] "DrugName"            "CellLineName"        "PrimaryTissue"
#  [4] "ExperimentNumber"    "TrtDuration"         "doublingtime"
#  [7] "log10Concentration"  "GR"                  "relative_cell_count"
# [10] "activities.mean"     "activities.std"      "Norm_CellLineName"
# [13] "Norm_DrugName"


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

raw_dt
message(paste0("Loading TREDataMapper"))
TREDataMapper <- CoreGx::TREDataMapper(rawdata=raw_dt)

CoreGx::rowDataMap(TREDataMapper) <- list(
    id_columns = (c("gCSI.treatmentid", "ExperimentNumber", "dose")),
    mapped_columns = c())

CoreGx::colDataMap(TREDataMapper) <- list(
    id_columns = c("gCSI.sampleid"),
    mapped_columns = c())

CoreGx::assayMap(TREDataMapper) <- list(
    raw = list(
        c("gCSI.treatmentid", "gCSI.sampleid", "ExperimentNumber", "dose"),
        c("viability")),
    summary = list(
        c("gCSI.treatmentid", "gCSI.sampleid", "ExperimentNumber", "dose"),
        c("GR", "activities.mean", "activities.std")))

gCSI_tre <- CoreGx::metaConstruct(TREDataMapper)
gCSI_tre





GRMetrics <- input_dt_list[["gCSI_GRmetrics_v1.3.tsv"]]
# names(GRMetrics)
# names(GRMetrics)
#  [1] "DrugName"              "CellLineName"          "PrimaryTissue"
#  [4] "ExperimentNumber"      "TrtDuration"           "doublingtime"
#  [7] "maxlog10Concentration" "GR_AOC"                "meanviability"
# [10] "GRmax"                 "Emax"                  "GRinf"
# [13] "GEC50"                 "GR50"                  "R_square_GR"
# [16] "pval_GR"               "flat_fit_GR"           "h_GR"
# [19] "N_conc"                "log10_conc_step"       "Norm_CellLineName"
# [22] "Norm_DrugName"         "GR_05uM_fit"


published_profiles <- unique(
    GRMetrics[,
    .(
    gCSI.treatmentid = Norm_DrugName,
    gCSI.sampleid = Norm_CellLineName,
    mean_viability = meanviability,
    GR_AOC,
    GRmax,
    Emax,
    GRinf,
    GEC50,
    GR50)])

CoreGx::assay(gCSI_tre, "profiles_published") <- published_profiles

qs::qsave(gCSI_tre, file = OUTPUT[[1]])

# tre_fit <- gCSI_tre[,1:6,] |> CoreGx::endoaggregate(
#     {  # the entire code block is evaluated for each group in our group by
#         # 1. fit a log logistic curve over the dose range
#         fit <- PharmacoGx::logLogisticRegression(dose, viability,
#             viability_as_pct=FALSE)
#         # 2. compute curve summary metrics
#         ic50 <- PharmacoGx::computeIC50(dose, Hill_fit=fit)
#         aac <- PharmacoGx::computeAUC(dose, Hill_fit=fit)
#         # 3. assemble the results into a list, each item will become a
#         #   column in the target assay.
#         list(
#             HS=fit[["HS"]],
#             E_inf = fit[["E_inf"]],
#             EC50 = fit[["EC50"]],
#             Rsq=as.numeric(unlist(attributes(fit))),
#             aac_recomputed=aac,
#             ic50_recomputed=ic50
#         )
#     },
#     assay="raw",
#     target="profiles_recomputed",
#     enlist=FALSE,  # this option enables the use of a code block for aggregation
#     by=c("gCSI.treatmentid", "gCSI.sampleid"),
#     nthread=6  # parallelize over multiple cores to speed up the computation
# )




