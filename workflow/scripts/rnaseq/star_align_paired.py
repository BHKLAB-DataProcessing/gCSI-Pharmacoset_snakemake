from snakemake.shell import shell

import tempfile

fq1 = snakemake.input.get("fq1", "")
fq2 = snakemake.input.get("fq2", "")
index = snakemake.input.get("index", "")
extra = snakemake.params.get("extra", "")

assert fq1 and fq2, "input-> fq1 and fq2 are required."

if fq1.endswith(".gz"):
    readcmd = "--readFilesCommand gunzip -c"
elif fq1.endswith(".bz2"):
    readcmd = "--readFilesCommand bunzip2 -c"
else:
    readcmd = ""

out_unmapped = snakemake.output.get("unmapped","")
if out_unmapped:
    out_unmapped = "--outReadsUnmapped Fastx"

if "--outSAMtype BAM SortedByCoordinate" in extra:
    stdout = "BAM_sortedByCoordinate > {snakemake.output.sorted_bam}"
elif "BAM_Unsorted" in extra:
    stdout = "BAM_Unsorted"
else:
    stdout = "SAM > {snakemake.output.algned_sam}"

log = snakemake.log_fmt_shell(stdout=False, stderr=True)

with tempfile.TemporaryDirectory() as tmpdir:
    shell(
        "STAR "
        " --runThreadN {snakemake.threads}"
        " --genomeDir {index}"
        " --readFilesIn {fq1} {fq2}"
        " {readcmd}"
        " {extra}"
        " {out_unmapped}"   
        " --outTmpDir {tmpdir}/STARtmp"
        " --outFileNamePrefix {tmpdir}/"
        " --outStd {stdout}"
        " {log}"
    )
    
    if snakemake.output.get("reads_per_gene"):
        shell("cat {tmpdir}/ReadsPerGene.out.tab > {snakemake.output.reads_per_gene:q}")
    if snakemake.output.get("chim_junc"):
        shell("cat {tmpdir}/Chimeric.out.junction > {snakemake.output.chim_junc:q}")
    if snakemake.output.get("sj"):
        shell("cat {tmpdir}/SJ.out.tab > {snakemake.output.sj:q}")
    if snakemake.output.get("log"):
        shell("cat {tmpdir}/Log.out >> {snakemake.output.log:q}")
    if snakemake.output.get("log_progress"):
        shell("cat {tmpdir}/Log.progress.out > {snakemake.output.log_progress:q}")
    if snakemake.output.get("log_final"):
        shell("cat {tmpdir}/Log.final.out > {snakemake.output.log_final:q}")
        
    unmapped = snakemake.output.get("unmapped")
    if unmapped:
        # SE
        if not fq2:
            unmapped = [unmapped]

        for i, out_unmapped in enumerate(unmapped, 1):
            cmd = "gzip -c" if out_unmapped.endswith("gz") else "cat"
            shell("{cmd} {tmpdir}/Unmapped.out.mate{i} > {out_unmapped}")
