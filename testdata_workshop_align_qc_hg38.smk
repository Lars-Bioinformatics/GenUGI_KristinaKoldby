# Run script from "/work/sduhumac/kristina/data/testdata_workshop/data"

from os import getcwd

SAMPLES, = glob_wildcards("{sample}_R1_001.fastq.gz")
FAMNAME = getcwd().rsplit("/",2)[1]

resource_path = "/work/sdukoldby/resources/hg38/"
ref = resource_path + "Homo_sapiens_assembly38.fasta"
dbsnp= resource_path + "Homo_sapiens_assembly38.dbsnp138.vcf"
mills_1000g=resource_path + "Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
phase1_1000g=resource_path + "1000G_phase1.snps.high_confidence.hg38.vcf.gz"
bed = resource_path + "workshop_chr21_regions_g1k_hg38.bed"
interval_list = resource_path + "workshop_chr21_regions_g1k_hg38_edited2.interval_list"

### Can be run without bed file

rule all:
  input:
    expand("{sample}.recalibrated.bam", sample=SAMPLES),
    expand("{sample}.recalibrated.bai", sample=SAMPLES),
    expand("quality_control/{famname}_multiqc_report.html", famname = FAMNAME)


rule fastqc:
    input:
        "{sample}_R1_001.fastq.gz",
        "{sample}_R2_001.fastq.gz",
    output:
        html_1 = "quality_control/fastqc/{sample}_R1_001_fastqc.html",
        zip_1 = "quality_control/fastqc/{sample}_R1_001_fastqc.zip",
        html_2 = "quality_control/fastqc/{sample}_R2_001_fastqc.html",
        zip_2 = "quality_control/fastqc/{sample}_R2_001_fastqc.zip",
    shell:
        """
        fastqc {input} -o "quality_control/fastqc/"
        """

rule bwa_mem:
    input:
        "{sample}_R1_001.fastq.gz",
        "{sample}_R2_001.fastq.gz"
    output:
        bam = temp("{sample}.sorted.bam"),
        bai = temp("{sample}.sorted.bai")
    threads: 24
    log:
        "{sample}.bwa_mem.log"
    params:
        rgid = "{sample}",
        rglb = "Illumina",
        rgsm = "{sample}",
        rgpl = "HiSeq",
        rgpu = "C2FAGACXX"
    shell:
        """
        bwa mem {ref} {input} \
        -R \"@RG\\tID:{params.rgid}\\tLB:{params.rglb}\\tSM:{params.rgsm}\\tPL:{params.rgpl}\\tPU:{params.rgpu}\" \
        -M \
        -t {threads} | \
        gatk --java-options -Xmx12G SortSam \
        --INPUT=/dev/stdin \
        --OUTPUT={output.bam} \
        --VALIDATION_STRINGENCY=LENIENT \
        --SORT_ORDER=coordinate \
        --CREATE_INDEX=TRUE
        """

rule markdup:
    input:
        bam = "{sample}.sorted.bam",
        bai = "{sample}.sorted.bai"
    output:
        bam = temp("{sample}.marked.bam"),
        bai = temp("{sample}.marked.bai"),
        metrics = "quality_control/{sample}_duplicate_metrics.txt"
    threads: 24
    log:
        "{sample}.markdup.log"
    shell:
        """
        gatk --java-options -Xmx12G MarkDuplicates \
        --INPUT={input.bam} \
        --OUTPUT={output.bam} \
        --METRICS_FILE={output.metrics} \
        --REMOVE_DUPLICATES=FALSE \
        --OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
        --CREATE_INDEX=TRUE
        """

# OPTICAL_DUPLICATE_PIXEL_DISTANCE:
# Should be set to 2500 for patterned flowcells like Illumina Hiseq 3000/4000 and Novaseq 6000
# and set to 100 for unpatterned flowcells like HiSeq 2500

rule BaseRecalibrator:
    input:
        bam = "{sample}.marked.bam",
        bai = "{sample}.marked.bai",
        bed = {bed}
    output:
        "quality_control/{sample}_pre_recalibration.grp"
    shell:
        """
        gatk --java-options -Xmx12G BaseRecalibrator \
        -R {ref} \
        -I {input.bam} \
        --known-sites {dbsnp} \
        --known-sites {mills_1000g} \
        --known-sites {phase1_1000g} \
        -O {output} \
        -L {input.bed}
        """

### Can be run without bed file


rule ApplyBQSR:
    input:
        bam = "{sample}.marked.bam",
        bai = "{sample}.marked.bai",
        recal = "quality_control/{sample}_pre_recalibration.grp",
        bed = {bed}
    output:
        bam = "{sample}.recalibrated.bam",
        bai = "{sample}.recalibrated.bai"
    shell:
        """
        gatk --java-options -Xmx12G ApplyBQSR \
        -R {ref} \
        -I {input.bam} \
        --bqsr {input.recal} \
        -O {output.bam} \
        -L {input.bed}
        """

rule recalibration_table:
    input:
        bam = "{sample}.recalibrated.bam",
        bai = "{sample}.recalibrated.bai",
        bed = {bed}
    output:
        "quality_control/{sample}_post_recalibration.grp"
    shell:
        """
        gatk --java-options -Xmx12G BaseRecalibrator \
        -R {ref} \
        -I {input.bam} \
        --known-sites {dbsnp} \
        --known-sites {mills_1000g} \
        --known-sites {phase1_1000g} \
        -O {output} \
        -L {input.bed}
        """


rule CollectHsMetrics:
    input:
        bam = "{sample}.recalibrated.bam",
        bai = "{sample}.recalibrated.bai",
        interval = {interval_list}
    output:
        "quality_control/{sample}_HsMetrics.txt"
    shell:
        """
        gatk --java-options -Xmx12G CollectHsMetrics \
        --INPUT={input.bam} \
        --REFERENCE_SEQUENCE={ref} \
        --OUTPUT={output} \
        --BAIT_INTERVALS={input.interval} \
        --TARGET_INTERVALS={input.interval}
        """

rule qualimap:
    input:
        "{sample}.recalibrated.bam"
    output:
        html = "quality_control/{sample}/qualimapReport.html"
    params:
        outdir = directory("quality_control/{sample}")
    shell:
        """
        qualimap bamqc \
        -bam {input} \
        -nt 24 \
        -c \
        -sd \
        -gff {bed} \
        -outdir {params.outdir} \
        --java-mem-size=4G
        """



rule multiqc:
    input:
        expand("quality_control/fastqc/{sample}_R1_001_fastqc.html", sample=SAMPLES),
        expand("quality_control/fastqc/{sample}_R1_001_fastqc.zip", sample=SAMPLES),
        expand("quality_control/{sample}_pre_recalibration.grp", sample=SAMPLES),
        expand("quality_control/{sample}/qualimapReport.html", sample=SAMPLES),
        expand("quality_control/{sample}_HsMetrics.txt", sample=SAMPLES),
        expand("quality_control/fastqc/{sample}_R2_001_fastqc.html", sample=SAMPLES),
        expand("quality_control/fastqc/{sample}_R2_001_fastqc.zip", sample=SAMPLES),
        expand("quality_control/{sample}_post_recalibration.grp", sample=SAMPLES),
        expand("quality_control/{sample}_duplicate_metrics.txt", sample=SAMPLES),
    output:
        html = "quality_control/{famname}_multiqc_report.html"
    params:
        indir = directory("quality_control"),
        config = resource_path + "multiqc_config.yaml"
    shell:
        """
        multiqc {params.indir} \
        -n {output.html} \
        -c {params.config}
        """
