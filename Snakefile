import pandas as pd
from os import path, listdir, makedirs
from os.path import join, exists, isfile, isdir
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider

HTTP = HTTPRemoteProvider()

configfile: "config/config.yaml"

REFERENCE_DIR = "references"

ref_build = ["GRCh38"]
ref_version = "45"
RSEM_VERSION = "1.3.0"
STAR_VERSION = "2.7.9a"
KALLISTO_version = "0.46.1"
SALMON_version = "1.8.0"
# metadata = "metadata"

# sample_df = pd.read_csv("metadata/rnaseq/sample_file.csv")

# SAMPLES = list(set(sample_df.sample_alias.to_list()))

include: "workflow/rules/downloadData.smk"
include: "workflow/rules/metadata.smk"
# include: "workflow/rules/rnaseq.smk" # processed on H4H compute server 

rule build_PharmacoSet:
    input:
        treatmentMetadata = "results/data/metadata/treatmentMetadata_annotated.tsv",
        sampleMetadata = "results/data/metadata/sampleMetadata_annotated.tsv",
        mae = "results/data/gCSI_MultiAssayExperiment.RDS",
        tre = "results/data/gCSI_TreatmentResponseExperiment.RDS"
    output:
        pharmacoSet = "results/data/gCSI_PharmacoSet.RDS"
    log:
        "logs/rnaseq/build_PharmacoSet.log"
    threads:
        8
    script:
        "workflow/scripts/build_PharmacoSet.R"


rule create_MultiAssayExperiment:
    input:
        rnaseq = HTTP.remote(config["molecularProfiles"]["rnaseq"]["url"]),
        sampleMetadata = "results/data/metadata/sampleMetadata_annotated.tsv"
    output:
        mae = "results/data/gCSI_MultiAssayExperiment.RDS"
    log:
        "logs/rnaseq/create_MultiAssayExperiment.log"
    threads:
        8
    script:
        "workflow/scripts/create_MultiAssayExperiment.R"


rule create_treatmentResponseExperiment:
    input:
        tsv_path = "rawdata/TreatmentResponse_and_Metadata/gCSI_GRdata_v1.3.tsv.tar.gz",
        treatmentMetadata = "results/data/metadata/treatmentMetadata_annotated.tsv",
        sampleMetadata = "results/data/metadata/sampleMetadata_annotated.tsv"
    output:
        tre = "results/data/gCSI_TreatmentResponseExperiment.RDS"
    log:
        "logs/rnaseq/create_treatmentResponseExperiment.log"
    threads:
        30
    script:
        "workflow/scripts/create_treatmentResponseExperiment.R"


