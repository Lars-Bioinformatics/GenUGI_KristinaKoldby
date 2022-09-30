WORK="/work/sdukoldby/data/G45-2016_genugi/190325_A00653_0016_AHGL2LDSXX/BaseCalls/hg38_withoutMinusM"
# WORK="/work/sdukoldby/data/G45-2016_genugi/190325_A00653_0016_AHGL2LDSXX/BaseCalls/hg38"
INPUT_BAM = WORK
PON_DIR = WORK+"/GRIDSS_joint_calling"
OUTPUT = WORK+"/GRIDSS"
GRIDSS_VERSION = "2.8.1"
SCRIPTS = WORK+"/GRIDSS_joint_calling/gridss-"+GRIDSS_VERSION

configfile: WORK+"/samples_full_names.yaml"

# input sample names
# SAMPLES, = glob_wildcards(INPUT_BAM+"/{sample}.recalibrated.bam")
# SAMPLES = sorted(SAMPLES)
SAMPLES = [config[n]["tumor"] for n in config] + [config[n]["normal"] for n in config]
print(SAMPLES)

NORMALS = [config[n]["normal"] for n in config]
print(NORMALS)

PON = "small_pon"

# output sample names
# PROJECT = "G45-ECV2"

# Resources
ref = "/work/sdukoldby/resources/hg38/hg38_new_names/Homo_sapiens_assembly38.fasta"
ref_r = "BSgenome.Hsapiens.UCSC.hg38"
res = "/work/sdukoldby/resources/hg38"
# mem = 12
mem = 31

onstart:
    shell("mkdir -p " + OUTPUT)

rule all:
    input:
        # [expand(OUTPUT+"/{tumor}_vs_{normal}__{pon}.filtered.vcf.gz", tumor=config[s]["tumor"], normal=config[s]["normal"], pon=PON) for s in config]
        # [expand(OUTPUT+"/{tumor}_vs_{normal}_unfiltered.vcf", tumor=config[s]["tumor"], normal=config[s]["normal"]) for s in config]
        [expand(OUTPUT+"/{tumor}_vs_{normal}__{pon}.somatic.annotated.vcf", tumor=config[s]["tumor"], normal=config[s]["normal"], pon=PON) for s in config]


rule gridss_pair:
    input:
        normal=INPUT_BAM+"/{normal}.recalibrated.bam",
        tumor=INPUT_BAM+"/{tumor}.recalibrated.bam",
        blacklist=res+"/ENCODE_DAC_blacklist/ENCFF419RSJ.bed"
    output:
        vcf=OUTPUT+"/{tumor}_vs_{normal}_unfiltered.vcf"
    threads: 8
    shell:
        """
        mkdir -p {OUTPUT}/{wildcards.tumor}_intermediate_results

        {SCRIPTS}/gridss.sh \
        -t {threads} \
        -b {input.blacklist} \
        -r {ref} \
        -w {OUTPUT}/{wildcards.tumor}_intermediate_results \
        -o {output.vcf} \
        -a {OUTPUT}/{wildcards.tumor}_vs_{wildcards.normal}_gridss_assembly.bam \
        -j {SCRIPTS}/gridss-{GRIDSS_VERSION}-gridss-jar-with-dependencies.jar \
        --jvmheap {mem}g \
        {input.normal} {input.tumor}
        """

rule gridss_normals:
    input:
        bam=expand(INPUT_BAM+"/{normal}.recalibrated.bam", normal=NORMALS),
        blacklist=res+"/ENCODE_DAC_blacklist/ENCFF419RSJ.bed"
    output:
        vcf=PON_DIR+"/{pon}dir/{pon}_gridss.vcf"
    threads: 8
    shell:
        """
        mkdir -p {OUTPUT}/{wildcards.pon}dir

        {SCRIPTS}/gridss.sh \
        -t {threads} \
        -b {input.blacklist} \
        -r {ref} \
        -w {PON_DIR}/intermediate_files \
        -o {output.vcf} \
        -a {PON_DIR}/intermediate_files/{wildcards.pon}.gridss.{PROJECT}.assembly.bam \
        -j {SCRIPTS}/gridss-{GRIDSS_VERSION}-gridss-jar-with-dependencies.jar \
        --jvmheap {mem}g \
        {input.bam}
        """

rule panel_of_normals:
    input:
        vcf=PON_DIR+"/{pon}dir/{pon}_gridss.vcf"
    output:
        bedpe=PON_DIR+"/{pon}dir/gridss_pon_breakpoint.bedpe",
        bed=PON_DIR+"/{pon}dir/gridss_pon_single_breakend.bed"
    shell:
        """
        java -Xmx{mem}g \
        	-cp {SCRIPTS}/gridss-{GRIDSS_VERSION}-gridss-jar-with-dependencies.jar \
        	gridss.GeneratePonBedpe \
        	INPUT={input.vcf} \
        	O={output.bedpe} \
        	SBO={output.bed} \
        	REFERENCE_SEQUENCE={ref}
        """

# Requires installed BsGenome e.g.
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
#
# BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
# BiocManager::install("VariantAnnotation")
# BiocManager::install("StructuralVariantAnnotation")
# Other R packages required: argparser, tidyverse, stringdist
rule filter_gridss:
    input:
        vcf=OUTPUT+"/{tumor}_vs_{normal}_unfiltered.vcf",
        bedpe=PON_DIR+"/{pon}dir/gridss_pon_breakpoint.bedpe",
        bed=PON_DIR+"/{pon}dir/gridss_pon_single_breakend.bed"
    output:
        vcf=OUTPUT+"/{tumor}_vs_{normal}__{pon}.somatic.vcf.gz"
    params:
        vcf=OUTPUT+"/{tumor}_vs_{normal}__{pon}.somatic.vcf"
    shell:
        """
        mkdir -p {OUTPUT}/plots

        Rscript {SCRIPTS}/gridss_somatic_filter.R \
        --pondir {PON_DIR}/{wildcards.pon}dir \
        --ref {ref_r} \
        --input {input.vcf} \
        --output {params.vcf} \
        --plotdir {OUTPUT}/plots \
        --scriptdir {SCRIPTS}
    
        rename bgz gz *bgz*
        """

rule annotate_svclass:
    input:
        vcf=OUTPUT+"/{tumor}_vs_{normal}__{pon}.somatic.vcf.gz"
    output:
        vcf=OUTPUT+"/{tumor}_vs_{normal}__{pon}.somatic.annotated.vcf",
        bed=OUTPUT+"/{tumor}_vs_{normal}__{pon}.somatic.annotated.bed"
    shell:
        """
        Rscript {SCRIPTS}/simple-event-annotation.R \
        --in_vcf {input.vcf} \
        --out_vcf {output.vcf} \
        --out_bed {output.bed}
        """


###############################################################################
### OLD STUFF
###############################################################################
rule gridss_conda:
    input:
        normal=INPUT_BAM+"/{normal}.recalibrated.bam",
        tumor=INPUT_BAM+"/{tumor}.recalibrated.bam",
        blacklist=res+"ENCODE_DAC_blacklist/ENCFF419RSJ.bed"
    output:
        vcf=OUTPUT+"/{tumor}_vs_{normal}_unfiltered_conda.vcf"
    threads: 8
    shell:
        """
        mkdir -p {OUTPUT}/{wildcards.tumor}_intermediate_results
        
        gridss gridss.CallVariants -Xmx{mem}g \
            WORKER_THREADS={threads} \
            REFERENCE_SEQUENCE={ref} \
            BLACKLIST={input.blacklist} \
            INPUT={input.normal} \
            INPUT={input.tumor} \
            OUTPUT={output.vcf} \
            ASSEMBLY={OUTPUT}/{wildcards.tumor}_vs_{wildcards.normal}_gridss_assembly.bam \
            WORKING_DIR={OUTPUT}/{wildcards.tumor}_intermediate_results
        """