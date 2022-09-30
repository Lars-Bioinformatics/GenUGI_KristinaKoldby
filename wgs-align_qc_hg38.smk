# Run script from "/work/sduhumac/kristina/data/testdata_workshop/data"

from os import getcwd
#SAMPLES = "G35-16ouh013390_illumina-truseq-genome_HLYWNBGXY"
SAMPLES, = glob_wildcards("{sample}_R1_001.fastq.gz")
FAMNAME = getcwd().rsplit("/",1)[1]

resource_path = "/work/sdukoldby/resources/hg38/"
ref = resource_path + "Homo_sapiens_assembly38.fasta"
dbsnp= resource_path + "Homo_sapiens_assembly38.dbsnp138.vcf"
mills_1000g=resource_path + "Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
phase1_1000g=resource_path + "1000G_phase1.snps.high_confidence.hg38.vcf.gz"


rule all:
  input:
    expand("hg38_withoutMinusM/{sample}.recalibrated.bam", sample=SAMPLES),
    expand("hg38_withoutMinusM/{sample}.recalibrated.bai", sample=SAMPLES),
    expand("hg38_withoutMinusM/quality_control/{famname}_multiqc_report.html", famname = FAMNAME)


rule fastqc:
    input:
        "{sample}_R1_001.fastq.gz",
        "{sample}_R2_001.fastq.gz",
    output:
        html_1 = "hg38_withoutMinusM/quality_control/fastqc/{sample}_R1_001_fastqc.html",
        zip_1 = "hg38_withoutMinusM/quality_control/fastqc/{sample}_R1_001_fastqc.zip",
        html_2 = "hg38_withoutMinusM/quality_control/fastqc/{sample}_R2_001_fastqc.html",
        zip_2 = "hg38_withoutMinusM/quality_control/fastqc/{sample}_R2_001_fastqc.zip",
    shell:
        """
        fastqc {input} -o "hg38_withoutMinusM/quality_control/fastqc/"
        """

rule bwa_mem:
    input:
        "{sampleid}_{protocol}_{flowcell}_R1_001.fastq.gz",
        "{sampleid}_{protocol}_{flowcell}_R2_001.fastq.gz"
    output:
        bam = temp("hg38_withoutMinusM/{sampleid}_{protocol}_{flowcell}.sorted.bam"),
        bai = temp("hg38_withoutMinusM/{sampleid}_{protocol}_{flowcell}.sorted.bai")
    threads: 24
    log:
        "hg38_withoutMinusM/{sampleid}_{protocol}_{flowcell}.bwa_mem.log"
    params:
        rgid = "{sampleid}_{protocol}_{flowcell}",
        rglb = "{protocol}",
        rgsm = "{sampleid}_{protocol}_{flowcell}",
        rgpl = "illumina",
        rgpu = "{flowcell}"
    shell:
        """
        bwa mem {ref} {input} \
        -R \"@RG\\tID:{params.rgid}\\tLB:{params.rglb}\\tSM:{params.rgsm}\\tPL:{params.rgpl}\\tPU:{params.rgpu}\" \
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
        bam = "hg38_withoutMinusM/{sample}.sorted.bam",
        bai = "hg38_withoutMinusM/{sample}.sorted.bai"
    output:
        bam = temp("hg38_withoutMinusM/{sample}.marked.bam"),
        bai = temp("hg38_withoutMinusM/{sample}.marked.bai"),
        metrics = "hg38_withoutMinusM/quality_control/{sample}_duplicate_metrics.txt"
    threads: 24
    log:
        "hg38_withoutMinusM/{sample}.markdup.log"
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
        bam = "hg38_withoutMinusM/{sample}.marked.bam",
        bai = "hg38_withoutMinusM/{sample}.marked.bai",
    output:
        "hg38_withoutMinusM/quality_control/{sample}_pre_recalibration.grp"
    shell:
        """
        gatk --java-options -Xmx12G BaseRecalibrator \
        -R {ref} \
        -I {input.bam} \
        --known-sites {dbsnp} \
        --known-sites {mills_1000g} \
        --known-sites {phase1_1000g} \
        -O {output} \
        """

### Can be run without bed file


rule ApplyBQSR:
    input:
        bam = "hg38_withoutMinusM/{sample}.marked.bam",
        bai = "hg38_withoutMinusM/{sample}.marked.bai",
        recal = "hg38_withoutMinusM/quality_control/{sample}_pre_recalibration.grp",
    output:
        bam = "hg38_withoutMinusM/{sample}.recalibrated.bam",
        bai = "hg38_withoutMinusM/{sample}.recalibrated.bai"
    shell:
        """
        gatk --java-options -Xmx12G ApplyBQSR \
        -R {ref} \
        -I {input.bam} \
        --bqsr {input.recal} \
        -O {output.bam} \
        """

rule recalibration_table:
    input:
        bam = "hg38_withoutMinusM/{sample}.recalibrated.bam",
        bai = "hg38_withoutMinusM/{sample}.recalibrated.bai",
    output:
        "hg38_withoutMinusM/quality_control/{sample}_post_recalibration.grp"
    shell:
        """
        gatk --java-options -Xmx12G BaseRecalibrator \
        -R {ref} \
        -I {input.bam} \
        --known-sites {dbsnp} \
        --known-sites {mills_1000g} \
        --known-sites {phase1_1000g} \
        -O {output} \
        """


rule CollectWgsMetrics:
    input:
        bam = "hg38_withoutMinusM/{sample}.recalibrated.bam",
        bai = "hg38_withoutMinusM/{sample}.recalibrated.bai",
    output:
        "hg38_withoutMinusM/quality_control/{sample}_WgsMetrics.txt"
    shell:
        """
        gatk --java-options -Xmx12G CollectWgsMetrics \
        --INPUT={input.bam} \
        --REFERENCE_SEQUENCE={ref} \
        --OUTPUT={output} \
        """

rule qualimap:
    input:
        "hg38_withoutMinusM/{sample}.recalibrated.bam"
    output:
        html = "hg38_withoutMinusM/quality_control/{sample}/qualimapReport.html"
    params:
        outdir = directory("hg38_withoutMinusM/quality_control/{sample}")
    shell:
        """
        qualimap bamqc \
        -bam {input} \
        -nt 24 \
        -c \
        -sd \
        -outdir {params.outdir} \
        --java-mem-size=20G
        """



rule multiqc:
    input:
        expand("hg38_withoutMinusM/quality_control/fastqc/{sample}_R1_001_fastqc.html", sample=SAMPLES),
        expand("hg38_withoutMinusM/quality_control/fastqc/{sample}_R1_001_fastqc.zip", sample=SAMPLES),
        expand("hg38_withoutMinusM/quality_control/{sample}_pre_recalibration.grp", sample=SAMPLES),
        expand("hg38_withoutMinusM/quality_control/{sample}/qualimapReport.html", sample=SAMPLES),
        expand("hg38_withoutMinusM/quality_control/{sample}_WgsMetrics.txt", sample=SAMPLES),
        expand("hg38_withoutMinusM/quality_control/fastqc/{sample}_R2_001_fastqc.html", sample=SAMPLES),
        expand("hg38_withoutMinusM/quality_control/fastqc/{sample}_R2_001_fastqc.zip", sample=SAMPLES),
        expand("hg38_withoutMinusM/quality_control/{sample}_post_recalibration.grp", sample=SAMPLES),
        expand("hg38_withoutMinusM/quality_control/{sample}_duplicate_metrics.txt", sample=SAMPLES),
    output:
        html = "hg38_withoutMinusM/quality_control/{famname}_multiqc_report.html"
    params:
        indir = directory("hg38_withoutMinusM/quality_control"),
        config = resource_path + "multiqc_config.yaml"
    shell:
        """
        multiqc {params.indir} \
        -n {output.html} \
        -c {params.config}
        """
