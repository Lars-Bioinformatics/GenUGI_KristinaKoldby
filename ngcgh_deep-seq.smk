##### NOT WORKING #####
### Snakemake-part of script seems to be okay, dry-run (-np) looks fine
### ngCGH produces an error when trying to run it in the snakemake-env
### Works fine in base env (without snakemake) - apparently ngCGH is incompatible with something in snakemake installation?
### Use script ngcgh_deep-seq.sh


## Run from folder with recalibrated bam files ##
configfile: "samples_vcf.yaml"
vcf_path = "/work/sdukoldby/data/G45-2016_genugi/exome_fastq/exome_fastq_merged_2nd_novaseq_run/hg38/Connor/mutect2_somatic_joint_calling/single_sample_vcfs/"

PAIR, = glob_wildcards(vcf_path + "{sample}.vcf.gz")
# PAIR = ["PC1-18-op-P10_tumor_tagseq-medexome-deep-seq"]
# PAIR, = glob_wildcards("{tumor}_tumor_tagseq-medexome.connor.recalibrated.bam")
# PAIR = ["ECV2-8-biopsi-I1"]
# PAIR = ["ECV2-29-biopsi-C1_tumor"]
# PAIR = [pair+"_tumor" for pair in PAIR]
print(PAIR)

onstart:
    shell("mkdir -p ngcgh")

# regions = "/work/sdukoldby/resources/hg38/MedExome_hg38_capture_targets.bed"
window_size = 1000

rule all:
    input:
        [expand("ngcgh/{tumor}_vs_{normal}_medexome-deep-seq_ngcgh_w{window_size}.txt",
                    window_size=window_size,
        			normal=config[pair]["normal"],
        			tumor=config[pair]["tumor"]) for pair in PAIR]

rule ngcgh_deep:
    input:
        normal = "{normal}_normal_tagseq-medexome.connor.recalibrated.bam",
        tumor = "{tumor}_tumor_tagseq-medexome-deep-seq.connor.recalibrated.bam",
    output:
        "ngcgh/{tumor}_vs_{normal}_medexome-deep-seq_ngcgh_w{window_size}.txt"
    shell:
        """
        ngCGH \
        -w {window_size} \
        -o {output} \
        {input.normal} {input.tumor}
        """

rule ngcgh_std:
    input:
        normal = "{normal}_normal_tagseq-medexome.connor.recalibrated.bam",
        tumor = "{tumor}_tumor_tagseq-medexome.connor.recalibrated.bam",
    output:
        "ngcgh/{tumor}_vs_{normal}_medexome-deep-seq_ngcgh_w{window_size}.txt"
    shell:
        """
        ngCGH \
        -w {window_size} \
        -o {output} \
        {input.normal} {input.tumor}
        """
