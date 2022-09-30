WORK="/work/sdukoldby/data/G45-2016_genugi/190325_A00653_0016_AHGL2LDSXX/BaseCalls/"
INPUT = WORK+"hg38/"
OUTPUT = WORK+"cgp_bam/"

SAMPLES, = glob_wildcards(INPUT+"/{sample}.recalibrated.bam")

# Resources - paths inside docker
ref_fai = "/work/sdukoldby/resources/hg38/Homo_sapiens_assembly38.fasta.fai"

onstart:
    shell("mkdir -p " + OUTPUT)

rule all:
    input:
        expand(OUTPUT+"{sample}/{sample}.bam.bas", sample=SAMPLES)
        # expand(OUTPUT+"{sample}/{sample}.bam", sample=SAMPLES)

# rule link_gatk_bam:
#     input:
#         bam=INPUT+"{sample}.recalibrated.bam",
#         bai=INPUT+"{sample}.recalibrated.bai"
#     output:
#         bam=OUTPUT+"{sample}/{sample}.bam",
#         bai=OUTPUT+"{sample}/{sample}.bai"
#     shell:
#         """
#         mkdir -p {OUTPUT}{wildcards.sample}
#         ln {input.bam} {output.bam}
#         ln {input.bai} {output.bai}
#         """

rule bam_stats:
    input:
        bam=OUTPUT+"{sample}/{sample}.bam"
    output:
        bas=OUTPUT+"{sample}/{sample}.bam.bas"
    shell:
        """
        python2 ~/udocker-1.1.4/udocker run \
        --user=pi \
        --volume=/work:/work \
        cgpwgs \
        bam_stats \
        -i {input} \
        -o {output} \
        -r {ref_fai} \
        -@ 12
        """

# snakemake -s /work/sdukoldby/scripts/cgpwgs_preCommands_lars.smk -j 999 --cluster "sbatch -A sdukoldby_slim --time 24:00:00"