# Run script from folder with fastq files.

from os import getcwd
import time

SAMPLES, = glob_wildcards("{sample}_R1_001.fastq.gz")
#SAMPLES = "PC1-10-blod_normal_tagseq-medexome"
FAMNAME = getcwd().rsplit("/",1)[1]

resource_path = "/work/sdukoldby/resources/hg38/"
ref = resource_path + "Homo_sapiens_assembly38.fasta"
dbsnp= resource_path + "Homo_sapiens_assembly38.dbsnp138.vcf"
mills_1000g=resource_path + "Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
phase1_1000g=resource_path + "1000G_phase1.snps.high_confidence.hg38.vcf.gz"
bed = resource_path + "MedExome_hg38_capture_targets.bed"
interval_list = resource_path + "MedExome_hg38_capture_targets.interval_list"

### Can be run without bed file

rule all:
  input:
    expand("hg38/MarkDuplicates/{sample}.recalibrated.bam", sample=SAMPLES),
    expand("hg38/MarkDuplicates/{sample}.recalibrated.bai", sample=SAMPLES),
    expand("hg38/MarkDuplicates/quality_control/{sample}.insert_size_histogram.pdf", sample=SAMPLES),
    expand("hg38/MarkDuplicates/quality_control/{sample}.insert_size_metrics.txt", sample=SAMPLES),
    expand("hg38/MarkDuplicates/quality_control/{famname}_multiqc_report.html", famname = FAMNAME),
    expand("hg38/MarkDuplicates/quality_control/{sample}_capture-targets_HsMetrics.txt", sample=SAMPLES)


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

rule bwa_mem:
#demultiplexing using bcl2fastq names fastq file like this:
#{sampleid}_{protocol}_{flowcell}_{samplenumber}_R1_001.fastq.gz
#where {samplenumber} is for 1 to n samples is S1 to Sn according to SampleSheet.csv
    input:
        "{sampleid}_{sampletype}_{protocol}_R1_001.fastq.gz",
        "{sampleid}_{sampletype}_{protocol}_R2_001.fastq.gz"
    output:
        bam = "hg38/{sampleid}_{sampletype}_{protocol}.sorted.bam",
        bai = "hg38/{sampleid}_{sampletype}_{protocol}.sorted.bai"
    threads: 24
    log:
        "hg38/{sampleid}_{sampletype}_{protocol}.bwa_mem.log"
    params:
        rgid = "{sampleid}_{sampletype}_{protocol}",
        rglb = "{protocol}",
        rgsm = "{sampleid}_{sampletype}_{protocol}",
        rgpl = "illumina",
        rgpu = "merged"
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
        bam = "hg38/{sample}.sorted.bam",
        bai = "hg38/{sample}.sorted.bai"
    output:
        bam = temp("hg38/MarkDuplicates/{sample}.marked.bam"),
        bai = temp("hg38/MarkDuplicates/{sample}.marked.bai"),
        metrics = "hg38/MarkDuplicates/quality_control/{sample}_markdup_metrics.txt"
    threads: 24
    log:
        "hg38/MarkDuplicates/quality_control/{sample}.markdup.log"
    shell:
        """
        gatk --java-options -Xmx20G MarkDuplicates \
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
        bam = "hg38/MarkDuplicates/{sample}.marked.bam",
        bai = "hg38/MarkDuplicates/{sample}.marked.bai",
    output:
        "hg38/MarkDuplicates/quality_control/{sample}_pre_recalibration.grp"
    shell:
        """
        gatk --java-options -Xmx20G BaseRecalibrator \
        -R {ref} \
        -I {input.bam} \
        --known-sites {dbsnp} \
        --known-sites {mills_1000g} \
        --known-sites {phase1_1000g} \
        -O {output}
        """
#        -L {input.bed}
#bed = {bed}
### Can be run without bed file


rule ApplyBQSR:
    input:
        bam = "hg38/MarkDuplicates/{sample}.marked.bam",
        bai = "hg38/MarkDuplicates/{sample}.marked.bai",
        recal = "hg38/MarkDuplicates/quality_control/{sample}_pre_recalibration.grp",
    output:
        bam = "hg38/MarkDuplicates/{sample}.recalibrated.bam",
        bai = "hg38/MarkDuplicates/{sample}.recalibrated.bai"
    shell:
        """
        gatk --java-options -Xmx20G ApplyBQSR \
        -R {ref} \
        -I {input.bam} \
        --bqsr {input.recal} \
        -O {output.bam}
        """
#-L {input.bed}
#bed = {bed}

rule recalibration_table:
    input:
        bam = "hg38/MarkDuplicates/{sample}.recalibrated.bam",
        bai = "hg38/MarkDuplicates/{sample}.recalibrated.bai",
    output:
        "hg38/MarkDuplicates/quality_control/{sample}_post_recalibration.grp"
    shell:
        """
        gatk --java-options -Xmx20G BaseRecalibrator \
        -R {ref} \
        -I {input.bam} \
        --known-sites {dbsnp} \
        --known-sites {mills_1000g} \
        --known-sites {phase1_1000g} \
        -O {output}
        """

#-L {input.bed}
#bed = {bed}

rule CollectHsMetrics:
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

rule CollectInsertSizeMetrics:
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


rule qualimap:
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



rule multiqc:
    input:
        expand("hg38/quality_control/fastqc/{sample}_R1_001_fastqc.html", sample=SAMPLES),
        expand("hg38/quality_control/fastqc/{sample}_R1_001_fastqc.zip", sample=SAMPLES),
        expand("hg38/MarkDuplicates/quality_control/{sample}_pre_recalibration.grp", sample=SAMPLES),
        expand("hg38/MarkDuplicates/quality_control/{sample}/qualimapReport.html", sample=SAMPLES),
        expand("hg38/MarkDuplicates/quality_control/{sample}_capture-targets_HsMetrics.txt", sample=SAMPLES),
        expand("hg38/quality_control/fastqc/{sample}_R2_001_fastqc.html", sample=SAMPLES),
        expand("hg38/quality_control/fastqc/{sample}_R2_001_fastqc.zip", sample=SAMPLES),
        expand("hg38/MarkDuplicates/quality_control/{sample}_post_recalibration.grp", sample=SAMPLES),
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

rule version_log:
    input:
        "hg38/MarkDuplicates/quality_control/{famname}_multiqc_report.html"
    output:
        log="hg38/logfiles/exom_tagseq_align_qc_hg38-version-log_%Y_%m_%d:%X.txt"
    shell:
        """
        echo 'Script: exom_tagseq_align_qc_hg38' > {output};
        conda list --export >> {output}
        """
