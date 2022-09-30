# Run script from folder with fastq files.

from os import getcwd
import time

SAMPLES, = glob_wildcards("{sample}_R1_001.fastq.gz")
#SAMPLES = "ECV2-4-plasma170503_tumor_tagseq-medexome"
#SAMPLES = "PC1-10-blod_normal_tagseq-medexome"
FAMNAME = getcwd().rsplit("/",1)[1]

resource_path = "/work/sdukoldby/resources/hg38/"
ref = resource_path + "Homo_sapiens_assembly38.fasta"
dbsnp= resource_path + "Homo_sapiens_assembly38.dbsnp138.vcf"
mills_1000g=resource_path + "Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
phase1_1000g=resource_path + "1000G_phase1.snps.high_confidence.hg38.vcf.gz"
bed = resource_path + "MedExome_hg38_capture_targets.bed"
interval_list = resource_path + "MedExome_hg38_capture_targets.interval_list"

rule all:
 input:
     expand("hg38/MarkDuplicates/quality_control/{sample}.insert_size_histogram.pdf", sample=SAMPLES),
     expand("hg38/MarkDuplicates/quality_control/{sample}.insert_size_metrics.txt", sample=SAMPLES),
     expand("hg38/Connor/quality_control/{sample}.insert_size_histogram.pdf", sample=SAMPLES),
     expand("hg38/Connor/quality_control/{sample}.insert_size_metrics.txt", sample=SAMPLES),
     "hg38/logfiles/exom_qc_tagseq-version-log_%Y_%m_%d:%X.txt"


rule fastqc:
    input:
        "{sample}_R1_001.fastq.gz",
        "{sample}_R2_001.fastq.gz",
    output:
        html_1 = "hg38/quality_control/fastqc/{sample}_R1_001_fastqc.html",
        zip_1 = "hg38/quality_control/fastqc/{sample}_R1_001_fastqc.zip",
        html_2 = "hg38/quality_control/fastqc/{sample}_R2_001_fastqc.html",
        zip_2 = "hg38/quality_control/fastqc/{sample}_R2_001_fastqc.zip",
    shell:
        """
        fastqc {input} -o "hg38/quality_control/fastqc/"
        """

#####################################
####       Markdup pipeline      ####
#####################################

rule CollectHsMetrics_markdup:
    input:
        bam = "hg38/MarkDuplicates/{sample}.recalibrated.bam",
        bai = "hg38/MarkDuplicates/{sample}.recalibrated.bai",
        interval = {interval_list}
    output:
        "hg38/MarkDuplicates/quality_control/{sample}_capture-targets_HsMetrics.txt"
    shell:
        """
        gatk --java-options -Xmx20G CollectHsMetrics \
        --INPUT={input.bam} \
        --REFERENCE_SEQUENCE={ref} \
        --OUTPUT={output} \
        --BAIT_INTERVALS={input.interval} \
        --TARGET_INTERVALS={input.interval}
        """

rule CollectInsertSizeMetrics_markdup:
    input:
        bam = "hg38/MarkDuplicates/{sample}.recalibrated.bam"
    output:
        hist = "hg38/MarkDuplicates/quality_control/{sample}.insert_size_histogram.pdf",
        metrics = "hg38/MarkDuplicates/quality_control/{sample}.insert_size_metrics.txt"
    shell:
        """
        gatk --java-options -Xmx20G CollectInsertSizeMetrics \
        -I {input.bam} \
        -H {output.hist} \
        -O {output.metrics}
        """


rule qualimap_markdup:
    input:
        "hg38/MarkDuplicates/{sample}.recalibrated.bam"
    output:
        html = "hg38/MarkDuplicates/quality_control/{sample}/qualimapReport.html"
    params:
        outdir = directory("hg38/MarkDuplicates/quality_control/{sample}")
    shell:
        """
        qualimap bamqc \
        -bam {input} \
        -nt 6 \
        -c \
        -sd \
        -gff {bed} \
        -outdir {params.outdir} \
        --java-mem-size=10G
        """

rule multiqc_markdup:
    input:
        expand("hg38/quality_control/fastqc/{sample}_R1_001_fastqc.html", sample=SAMPLES),
        expand("hg38/quality_control/fastqc/{sample}_R1_001_fastqc.zip", sample=SAMPLES),
        expand("hg38/MarkDuplicates/quality_control/{sample}/qualimapReport.html", sample=SAMPLES),
        expand("hg38/MarkDuplicates/quality_control/{sample}_capture-targets_HsMetrics.txt", sample=SAMPLES),
        expand("hg38/quality_control/fastqc/{sample}_R2_001_fastqc.html", sample=SAMPLES),
        expand("hg38/quality_control/fastqc/{sample}_R2_001_fastqc.zip", sample=SAMPLES),
        expand("hg38/MarkDuplicates/quality_control/{sample}_markdup_metrics.txt", sample=SAMPLES)
    output:
        html = "hg38/MarkDuplicates/quality_control/{famname}_multiqc_report.html"
    params:
        indir1 = directory("hg38/MarkDuplicates/"),
        indir2 = directory("hg38/quality_control/fastqc/"),
        config = resource_path + "multiqc_config.yaml"
    shell:
        """
        multiqc {params.indir1} {params.indir2} \
        -n {output.html} \
        -c {params.config}
        """

#####################################
####       Connor pipeline      ####
#####################################

rule CollectHsMetrics_connor:
    input:
        bam = "hg38/Connor/{sample}.connor.recalibrated.bam",
        bai = "hg38/Connor/{sample}.connor.recalibrated.bai",
        interval = {interval_list}
    output:
        "hg38/Connor/quality_control/{sample}_capture-targets_HsMetrics.txt"
    shell:
        """
        gatk --java-options -Xmx20G CollectHsMetrics \
        --INPUT={input.bam} \
        --REFERENCE_SEQUENCE={ref} \
        --OUTPUT={output} \
        --BAIT_INTERVALS={input.interval} \
        --TARGET_INTERVALS={input.interval}
        """

rule CollectInsertSizeMetrics_connor:
    input:
        bam = "hg38/Connor/{sample}.connor.recalibrated.bam"
    output:
        hist = "hg38/Connor/quality_control/{sample}.insert_size_histogram.pdf",
        metrics = "hg38/Connor/quality_control/{sample}.insert_size_metrics.txt"
    shell:
        """
        gatk --java-options -Xmx20G CollectInsertSizeMetrics \
        -I {input.bam} \
        -H {output.hist} \
        -O {output.metrics}
        """


rule qualimap_connor:
    input:
        "hg38/Connor/{sample}.connor.recalibrated.bam"
    output:
        html = "hg38/Connor/quality_control/{sample}/qualimapReport.html"
    params:
        outdir = directory("hg38/Connor/quality_control/{sample}")
    shell:
        """
        qualimap bamqc \
        -bam {input} \
        -nt 24 \
        -c \
        -sd \
        -gff {bed} \
        -outdir {params.outdir} \
        --java-mem-size=40G
        """



rule multiqc_connor:
    input:
        expand("hg38/quality_control/fastqc/{sample}_R1_001_fastqc.html", sample=SAMPLES),
        expand("hg38/quality_control/fastqc/{sample}_R1_001_fastqc.zip", sample=SAMPLES),
        expand("hg38/Connor/quality_control/{sample}/qualimapReport.html", sample=SAMPLES),
        expand("hg38/Connor/quality_control/{sample}_capture-targets_HsMetrics.txt", sample=SAMPLES),
        expand("hg38/quality_control/fastqc/{sample}_R2_001_fastqc.html", sample=SAMPLES),
        expand("hg38/quality_control/fastqc/{sample}_R2_001_fastqc.zip", sample=SAMPLES)
    output:
        html = "hg38/Connor/quality_control/{famname}_connor_multiqc_report.html"
    params:
        indir1 = directory("hg38/Connor/"),
        indir2 = directory("hg38/quality_control/fastqc/"),
        config = resource_path + "multiqc_config.yaml"
    shell:
        """
        multiqc {params.indir1} {params.indir2} \
        -n {output.html} \
        -c {params.config}
        """

rule version_log:
    input:
        expand("hg38/Connor/quality_control/{famname}_connor_multiqc_report.html", famname=FAMNAME),
        expand("hg38/MarkDuplicates/quality_control/{famname}_multiqc_report.html", famname=FAMNAME)
    output:
        log="hg38/logfiles/exom_qc_tagseq-version-log_%Y_%m_%d:%X.txt"
    shell:
        """
        echo 'Script: exom_qc_tagseq' > {output};
        conda list --export >> {output}
        """
