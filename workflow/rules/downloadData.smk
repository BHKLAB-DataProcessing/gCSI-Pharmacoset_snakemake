from os.path import join
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider

HTTP = HTTPRemoteProvider()


rule downloadTreatmentResponseANDMetadata:
    input:
        metadata_tsv=HTTP.remote(config["treatmentResponse_and_Metadata"]["TSV_files"]),
        metadata_rds=HTTP.remote(config["treatmentResponse_and_Metadata"]["RDS_files"]),
    output:
        metadata_tsv="rawdata/TreatmentResponse_and_Metadata/gCSI_GRdata_v1.3.tsv.tar.gz",
        metadata_rds="rawdata/TreatmentResponse_and_Metadata/gCSI_GRdata_v1.3.rds.tar.gz",
    shell:
        """
        mv {input.metadata_tsv} {output.metadata_tsv} && \
        mv {input.metadata_rds} {output.metadata_rds}
        """


rule processTreatmentResponse:
    input:
        metadata_tsv="rawdata/TreatmentResponse_and_Metadata/gCSI_GRdata_v1.3.tsv.tar.gz",
        metadata_rds="rawdata/TreatmentResponse_and_Metadata/gCSI_GRdata_v1.3.rds.tar.gz",
    output:
        TRPreprocessed="procdata/TreatmentResponsePreprocessed.qs",
    script:
        "../scripts/treatmentResponse/processTreatmentResponse.R"


# rule downloadExpressionMetadata:
#     input:
#         HTTP.remote(config["molecularProfiles"]["expression"]["metadata"]["url"])
#     output:
#         join(metadata,"EGAD000010000725/EGAD000010000725_map.txt")
#     shell:
#         "wget -O {output} {input}"
