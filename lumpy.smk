


ref = "/work/sdukoldby/resources/hg38/Homo_sapiens_assembly38.fasta"
exclude = "/work/sdukoldby/resources/hg38/exclude.cnvnator_100bp.GRCh38.20170403.bed"

configfile: "samples_full_names.yaml"
outdir = "lumpy"


rule all:
    input:
        [expand(outdir+"/{tumor}_vs_{normal}-smoove.genotyped.vcf.gz",
            normal=config[s]["normal"],
            tumor=config[s]["tumor"]) for s in config]


rule lumpy_smoove:
    input:
        normal = "{normal}.recalibrated.bam",
        tumor = "{tumor}.recalibrated.bam"
    output:
        vcf = "lumpy/{tumor}_vs_{normal}-smoove.genotyped.vcf.gz"
    threads: 24
    shell:
        """
        smoove call \
        -x \
        --name {wildcards.tumor}_vs_{wildcards.normal} \
        --outdir {outdir}
        --exclude {exclude} \
        --fasta {ref} \
        -p {threads} \
        --genotype {input}
        """
