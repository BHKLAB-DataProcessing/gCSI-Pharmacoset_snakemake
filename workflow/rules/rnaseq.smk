
configfile: "workflow/config/config.yaml"


tool_dict = {
    "kallisto": "abundance.h5",
    "salmon": "quant.sf",
    "rsem": ["genes.results", "isoforms.results"]
}


rule create_rse_list_RSEM:
    input:
        quant_files = lambda wc:
            expand(
                "procdata/rnaseq/rsem_v{RSEM_VERSION}_STAR_v{STAR_VERSION}_{ref_build}.{ref_version}/{sample}/{sample}.{file}",
                RSEM_VERSION= wc.tool_version,
                STAR_VERSION= wc.tool2_version,
                ref_build=wc.ref_build,
                ref_version=wc.ref_version,
                sample=SAMPLES,  
                file = tool_dict["rsem"],
            ),
        annotation_file = lambda wc: expand(
            "{REFERENCE_DIR}/human/{ref_build}_v{ref_version}/annotation.gtf",
            REFERENCE_DIR=REFERENCE_DIR,
            ref_build=wc.ref_build, 
            ref_version=wc.ref_version, 
        )
    output:
        rse_list = "results/data/rnaseq/{tool}-{tool2}_v{tool_version}-v{tool2_version}/{ref_build}.{ref_version}/rse_list.RData",
    params:
        tool = "rsem",
        tool_version = RSEM_VERSION,
    log:
        "logs/rnaseq/{tool}-{tool2}_v{tool_version}-v{tool2_version}/{ref_build}.{ref_version}/create_rse_list.log"
    envmodules: "R/4.2.1"
    threads: 1
    script:
        "scripts/rnaseq/create_rse_list.R"

rule create_rse_list_SALMON_KALLISTO:
    input:
        quant_files = lambda wc:
            expand(
                "procdata/rnaseq/{TOOL}_v{TOOL_version}_{ref_build}.{ref_version}/{sample}/{file}",
                TOOL = wc.TOOL, 
                TOOL_version = wc.TOOL_version,
                ref_build = wc.ref_build, 
                ref_version = wc.ref_version,
                sample = SAMPLES,  
                file = tool_dict[wc.TOOL.strip()],
            ),
        annotation_file = lambda wc: expand(
            "{REFERENCE_DIR}/human/{ref_build}_v{ref_version}/annotation.gtf",
            REFERENCE_DIR=REFERENCE_DIR,
            ref_build=wc.ref_build, 
            ref_version=wc.ref_version, 
        )
    output:
        rse_list = "results/data/rnaseq/{TOOL}_v{TOOL_version}_{ref_build}.{ref_version}/rse_list.RData",
    params:
        tool = lambda wc: "kallisto" if "kallisto" in wc.TOOL else "salmon",
        tool_version = lambda wc: KALLISTO_version if "kallisto" in wc.TOOL else SALMON_version,
    log:
        "logs/rnaseq/{TOOL}_ALL_quant/v{TOOL_version}_{ref_build}.{ref_version}/{TOOL}_ALL_quant.txt"
    envmodules: "R/4.2.1"
    threads: 1
    script:
        "scripts/rnaseq/create_rse_list.R"


rule kallisto_quant:
    input:
        fastq=lambda wc: [
            "rawdata/rnaseq/fastq/{}_1_{}.rnaseq.fastq.gz".format(wc.sample, pair)
            for pair in [1, 2]
        ],
        index=lambda wc: "{}/human/{}_v{}/KALLISTO.v{}/kallisto-transcriptome.idx".format(
            REFERENCE_DIR, wc.ref_build, wc.ref_version, wc.KALLISTO_version
        ),
    output:
        files=multiext(
            "procdata/rnaseq/kallisto_v{KALLISTO_version}_{ref_build}.{ref_version}/{sample}/",
            "abundance.h5",
            "abundance.tsv",
            "run_info.json",
        ),
    log:
        "logs/rnaseq/kallisto_v{KALLISTO_version}_quant/{sample}_{ref_build}.{ref_version}_kallisto-quant.log",
    envmodules: f"kallisto/{KALLISTO_version}"
    threads: 8
    shell:
        """
        kallisto version 
        OUTPUT_DIR=$(dirname {output.files[1]}) 
        echo $OUTPUT_DIR 2>&1 
        mkdir -p $OUTPUT_DIR 2>&1 


        kallisto quant \
            --threads {threads} \
            --index {input.index} \
            --output-dir $OUTPUT_DIR \
            {input.fastq[0]} {input.fastq[1]} 2>&1 | tee {log};
        """

rule salmon_quant:
    input:
        fastq=lambda wc: [
            "rawdata/rnaseq/fastq/{}_1_{}.rnaseq.fastq.gz".format(wc.sample, pair)
            for pair in [1, 2]
        ],
        index=lambda wc: "{}/human/{}_v{}/SALMON.v{}_index".format(
            REFERENCE_DIR, wc.ref_build, wc.ref_version, wc.SALMON_version
        ),
        transcriptome=lambda wc: "{}/human/{}_v{}/transcriptome.fa".format(
            REFERENCE_DIR, wc.ref_build, wc.ref_version
        ),
    output:
        quant = "procdata/rnaseq/salmon_v{SALMON_version}_{ref_build}.{ref_version}/{sample}/quant.sf",
        lib = "procdata/rnaseq/salmon_v{SALMON_version}_{ref_build}.{ref_version}/{sample}/lib_format_counts.json",
        cmd_info = "procdata/rnaseq/salmon_v{SALMON_version}_{ref_build}.{ref_version}/{sample}/cmd_info.json",
    params:
        libType="A",
    log:
        "logs/rnaseq/salmon_v{SALMON_version}_quant/{sample}_{ref_build}.{ref_version}_salmon-quant.log",
    envmodules: f"salmon/{SALMON_version}"
    threads: 8
    shell:
        """
        salmon version 
        OUTPUT_DIR=$(dirname {output.quant}) 
        echo $OUTPUT_DIR 2>&1 
        mkdir -p $OUTPUT_DIR 2>&1 

        salmon quant \
            --threads {threads} \
            --index {input.index} \
            --libType {params.libType} \
            --output $OUTPUT_DIR \
            --validateMappings \
            --geneMap {input.transcriptome} \
            -1 {input.fastq[0]} \
            -2 {input.fastq[1]} 2>&1 | tee {log};
        """



rule star_align_paired:
    input:
        fq1="rawdata/rnaseq/fastq/{sample}_1_1.rnaseq.fastq.gz",
        fq2="rawdata/rnaseq/fastq/{sample}_1_2.rnaseq.fastq.gz",
        index=lambda wc: "{}/human/{}_v{}/STAR.v{}_index".format(
            REFERENCE_DIR, wc.ref_build, wc.ref_version, wc.STAR_VERSION
        ),
    output:
        # note there are other possible outputs, but I am not using them for now
        # algned_sam = "procdata/rnaseq/star_v{STAR_VERSION}_{ref_build}.{ref_version}/{sample}_aligned.sam",
        # sorted_bam="procdata/rnaseq/star_v{STAR_VERSION}_{ref_build}.{ref_version}/{sample}/Aligned.sortedByCoord.out.bam",
        aligned_same = temp("procdata/rnaseq/star_v{STAR_VERSION}_{ref_build}.{ref_version}/{sample}/Aligned.out.sam"),
        transcript_bam= temp("procdata/rnaseq/star_v{STAR_VERSION}_{ref_build}.{ref_version}/{sample}/Aligned.toTranscriptome.out.bam"),
        sj=temp("procdata/rnaseq/star_v{STAR_VERSION}_{ref_build}.{ref_version}/{sample}/SJ.out.tab"),
        logout="procdata/rnaseq/star_v{STAR_VERSION}_{ref_build}.{ref_version}/{sample}/Log.out",
        logfinal="procdata/rnaseq/star_v{STAR_VERSION}_{ref_build}.{ref_version}/{sample}/Log.final.out",
        logprogress="procdata/rnaseq/star_v{STAR_VERSION}_{ref_build}.{ref_version}/{sample}/Log.progress.out",
    log:
        "logs/rnaseq/star_v{STAR_VERSION}/align_paired/{ref_build}.{ref_version}/{sample}_aligned.log",
    threads: 16
    envmodules:f"STAR/{STAR_VERSION}",
    shell:
        """
        STAR \
            --runThreadN {threads} \
            --genomeDir {input.index} \
            --readFilesIn {input.fq1} {input.fq2} \
            --readFilesCommand gunzip -c \
            --outFileNamePrefix $(dirname {output.transcript_bam})/ \
            --outFilterMultimapNmax 1 \
            --outFilterMismatchNmax 3 \
            --outFilterMismatchNoverLmax 0.3 \
            --quantMode TranscriptomeSAM \
            --sjdbOverhang 100 \
            --alignIntronMax 500000 \
            --alignMatesGapMax 500000 2>&1 | tee {log};
        """

rule rsem_calculateExpression:
    input:
        bam=lambda wc: f"procdata/rnaseq/star_v{wc.STAR_VERSION}_{wc.ref_build}.{wc.ref_version}/{wc.sample}/Aligned.toTranscriptome.out.bam",
        reference = lambda wc: f"references/human/{wc.ref_build}_v{wc.ref_version}/RSEM.v{wc.RSEM_VERSION}_index/reference.seq",
    output:
        # this file contains per-gene quantification data for the sample
        genes_results="procdata/rnaseq/rsem_v{RSEM_VERSION}_STAR_v{STAR_VERSION}_{ref_build}.{ref_version}/{sample}/{sample}.genes.results",
        # this file contains per-transcript quantification data for the sample
        isoforms_results="procdata/rnaseq/rsem_v{RSEM_VERSION}_STAR_v{STAR_VERSION}_{ref_build}.{ref_version}/{sample}/{sample}.isoforms.results",
    params:
        extra="--seed 42 --paired-end --no-bam-output --time",
    threads: 16
    envmodules: 
        "perl",
        f"rsem/{RSEM_VERSION}"
    log:
        "logs/rnaseq/rsem_v{RSEM_VERSION}_STAR_v{STAR_VERSION}_{ref_build}.{ref_version}/{sample}_calculateExpression.log",
    shell:
        """
        rsem-calculate-expression \
            --num-threads {threads} \
            {params.extra} \
            --alignments {input.bam} \
            $(dirname {input.reference})/reference \
            $(dirname {output.genes_results})/{wildcards.sample} 2>&1 | tee {log};
        """

