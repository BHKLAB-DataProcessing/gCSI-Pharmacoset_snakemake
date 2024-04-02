configfile: "workflow/config/config.yaml"

AnnotationGxDocker = config["containers"]["annotationGx"]

rule annotate_treatmentMetadata:
    input:
        rawTreatmentMetadata = "procdata/metadata/treatmentMetadata.tsv",
    output:
        treatmentMetadata = "results/data/metadata/treatmentMetadata_annotated.tsv"
    log:
        "logs/metadata/treatmentMetadata_annotated.log"
    threads:
        8
    container:
        AnnotationGxDocker
    script:
        "workflow/scripts/metadata/annotate_treatmentMetadata.R"


rule annotate_sampleMetadata:
    input:
        rawSampleMetadata = "procdata/metadata/sampleMetadata.tsv",
        EGA_sampleMetadata = "metadata/rnaseq/E-MTAB-2706.sdrf.txt"
    output:
        sampleMetadata = "results/data/metadata/sampleMetadata_annotated.tsv"
    log:
        "logs/metadata/sampleMetadata_annotated.log"
    threads:
        8
    container:
        AnnotationGxDocker
    script:
        "workflow/scripts/metadata/annotate_sampleMetadata.R"

rule preprocess_metadata:
    input:
        rawData = "rawdata/TreatmentResponse_and_Metadata/gCSI_GRdata_v1.3.tsv.tar.gz"
    output:
        rawSampleMetadata = "procdata/metadata/sampleMetadata.tsv",
        rawTreatmentMetadata = "procdata/metadata/treatmentMetadata.tsv",
    log:
        "logs/metadata/preprocess_metadata_annotated.log"
    threads:
        1
    script:
        "workflow/scripts/metadata/preprocess_metadata.R"

