# Run script from folder with fastq files.

from os import getcwd

SAMPLES, = glob_wildcards("{sample}_R1_001.fastq.gz")
#SAMPLES = "PC1-10-blod_S29_takara-tagseq_BHGTHHDSXX"
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


rule Mutect2
    input:
        normal =
        tumor =
    output:

    shell:
