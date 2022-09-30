from os import getcwd

SAMPLES, = glob_wildcards("{sample}_R1_001.fastq.gz")
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
    expand("hg38/{sample}_capture-file_HsMetrics.txt", sample=SAMPLES),


rule CollectHsMetrics:
    input:
        bam = "hg38/{sample}.recalibrated.bam",
        bai = "hg38/{sample}.recalibrated.bai",
        interval = {interval_list}
    output:
        "hg38/quality_control/{sample}_capture-file_HsMetrics.txt"
    shell:
        """
        gatk --java-options -Xmx12G CollectHsMetrics \
        --INPUT={input.bam} \
        --REFERENCE_SEQUENCE={ref} \
        --OUTPUT={output} \
        --BAIT_INTERVALS={input.interval} \
        --TARGET_INTERVALS={input.interval}
        """
