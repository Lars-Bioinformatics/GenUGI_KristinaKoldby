WORK="/work/sdukoldby/data/G45-2016_genugi/190325_A00653_0016_AHGL2LDSXX/BaseCalls/hg38_withoutMinusM"
INPUT_BAM = WORK
NORMAL_BAM = WORK+"/GRIDSS_joint_calling/normal_bam"
OUTPUT = WORK+"/GRIDSS_joint_calling"
GRIDSS_VERSION = "2.8.1"
SCRIPTS = WORK+"/GRIDSS_joint_calling/gridss-"+GRIDSS_VERSION

configfile: OUTPUT+"/samples_full_names.yaml"

# input sample names - ordered by sample type i.e. tumor or normal
SAMPLES = [config[n]["tumor"] for n in config if n.startswith("G45")] + [config[n]["normal"] for n in config if n.startswith("G45")]
# input sample names - ordered by sample name
# SAMPLES = sorted([config[n]["tumor"] for n in config] + [config[n]["normal"] for n in config])
print(SAMPLES)

PONS = "small_pon"
# PONS = ["large_pon"]
# PONS = ["small_pon","large_pon"]
# NORMALS = [config[n]["normal"] for n in config]
# NORMALS = config[pon]
# print(NORMALS)

# output sample names
PROJECT = "G45-ECV2"

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
        expand(OUTPUT+"/{project}_joint_calling.filtered_{pon}.vcf.gz", project=PROJECT, pon=PONS)
        # expand(OUTPUT+"/{pon}dir/gridss_pon_breakpoint.bedpe", pon=PONS),
        # expand(OUTPUT+"/{pon}dir/gridss_pon_single_breakend.bed", pon=PON),
        # expand(OUTPUT+"/{project}_joint_calling.filtered.vcf", project=PROJECT),
        # expand(OUTPUT+"/{pon}_{project}.vcf", project=PROJECT, pon=pon)
        # [expand(OUTPUT+"/{tumor}_vs_{normal}_unfiltered.vcf", tumor=config[s]["tumor"], normal=config[s]["normal"]) for s in config]


rule gridss_preprocessing:
    input:
        bam=INPUT_BAM+"/{sample}.recalibrated.bam",
        blacklist=res+"/ENCODE_DAC_blacklist/ENCFF419RSJ.bed"
    output:
        bam=OUTPUT+"/intermediate_files/{sample}.recalibrated.bam.gridss.working/{sample}.recalibrated.bam.sv.bam"
    threads: 8
    shell:
        """
        mkdir -p {OUTPUT}/intermediate_files
        ulimit -n 4096

        {SCRIPTS}/gridss.sh \
        -t {threads} \
        -b {input.blacklist} \
        -r {ref} \
        -w {OUTPUT}/intermediate_files \
        -o {OUTPUT}/{wildcards.sample}.tmp.vcf \
        -a {OUTPUT}/intermediate_files/{wildcards.sample}.gridss.{PROJECT}.assembly.bam \
        -j {SCRIPTS}/gridss-{GRIDSS_VERSION}-gridss-jar-with-dependencies.jar \
        --jvmheap {mem}g \
        --steps PreProcess \
        {input.bam}
        """
        # export _JAVA_OPTIONS="-Dgridss.defensiveGC=true"

# Example:
# gridss.sh  -t 4 -b wgEncodeDacMapabilityConsensusExcludable.bed -r ../hg19.fa -w out -o out/gridss.full.chr12.1527326.DEL1024.vcf -a out/gridss.full.chr12.1527326.DEL1024.assembly.bam -j ../target/gridss-2.8.0-gridss-jar-with-dependencies.jar --jvmheap 8g chr12.1527326.DEL1024.bam
rule gridss:
    input:
        bam=expand(INPUT_BAM+"/{sample}.recalibrated.bam", sample=SAMPLES),
        bam_sv=expand(OUTPUT+"/intermediate_files/{sample}.recalibrated.bam.gridss.working/{sample}.recalibrated.bam.sv.bam", sample=SAMPLES),
        blacklist=res+"/ENCODE_DAC_blacklist/ENCFF419RSJ.bed"
    output:
        vcf=OUTPUT+"/{project}_joint_calling.unfiltered.vcf"
    threads: 8
    shell:
        """
        mkdir -p {OUTPUT}/intermediate_files
        ulimit -n 4096
        export _JAVA_OPTIONS="-Dgridss.gridss.output_to_temp_file=true"

        {SCRIPTS}/gridss.sh \
        -t {threads} \
        -b {input.blacklist} \
        -r {ref} \
        -w {OUTPUT}/intermediate_files \
        -o {output.vcf} \
        -a {OUTPUT}/intermediate_files/gridss.{PROJECT}.assembly.bam \
        -j {SCRIPTS}/gridss-{GRIDSS_VERSION}-gridss-jar-with-dependencies.jar \
        --jvmheap {mem}g \
        {input.bam}
        """

rule gridss_normals:
    input:
        # bam=expand(INPUT_BAM+"/{normal}.recalibrated.bam", normal=NORMALS),
        # bam=expand("{normal}", normal=NORMALS),
        # bam=lambda wildcards: expand("{normal}", normal=config[wildcards.pon]),
        # bam_sv=lambda wildcards: expand(OUTPUT+"/{{pon}}_intermediate_results/{normal}.bam.gridss.working/{normal}.bam.sv.bam", normal=config[wildcards.pon]),
        bam=lambda wildcards: expand(NORMAL_BAM+"/{normal}.recalibrated.bam", normal=config[wildcards.pon]),
        bam_sv=lambda wildcards: expand(OUTPUT+"/preprocessed_samples/{normal}.recalibrated.bam.gridss.working/{normal}.recalibrated.bam.sv.bam", normal=config[wildcards.pon]),
        blacklist=res+"/ENCODE_DAC_blacklist/ENCFF419RSJ.bed"
    output:
        vcf=OUTPUT+"/{pon}dir/{pon}_gridss.vcf"
    threads: 8
    shell:
        """
        mkdir -p {OUTPUT}/{wildcards.pon}dir
        ulimit -n 4096

        {SCRIPTS}/gridss.sh \
        -t {threads} \
        -b {input.blacklist} \
        -r {ref} \
        -w {OUTPUT}/intermediate_files \
        -o {output.vcf} \
        -a {OUTPUT}/intermediate_files/{wildcards.pon}.gridss.{PROJECT}.assembly.bam \
        -j {SCRIPTS}/gridss-{GRIDSS_VERSION}-gridss-jar-with-dependencies.jar \
        --jvmheap {mem}g \
        --steps Assemble,Call \
        {input.bam}
        """
        # export _JAVA_OPTIONS="-Dgridss.defensiveGC=true"

rule panel_of_normals:
    input:
        vcf=OUTPUT+"/{pon}dir/{pon}_gridss.vcf"
    output:
        bedpe=OUTPUT+"/{pon}dir/gridss_pon_breakpoint.bedpe",
        bed=OUTPUT+"/{pon}dir/gridss_pon_single_breakend.bed"
    shell:
        """
        java -Xmx8g \
        	-cp {SCRIPTS}/gridss-{GRIDSS_VERSION}-gridss-jar-with-dependencies.jar \
        	gridss.GeneratePonBedpe \
        	INPUT={input.vcf} \
        	O={output.bedpe} \
        	SBO={output.bed} \
        	REFERENCE_SEQUENCE={ref}
        """
        # Include filtering here - seen with at least 3 occurrences


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
        vcf=OUTPUT+"/{project}_joint_calling.unfiltered.vcf",
        bedpe=OUTPUT+"/{pon}dir/gridss_pon_breakpoint.bedpe",
        bed=OUTPUT+"/{pon}dir/gridss_pon_single_breakend.bed"
    output:
        vcf=OUTPUT+"/{project}_joint_calling.filtered_{pon}.vcf.gz"
    params:
        vcf=OUTPUT+"/{project}_joint_calling.filtered_{pon}.vcf"
    shell:
        """
        mkdir -p {OUTPUT}/plots

        Rscript {SCRIPTS}/gridss_somatic_filter.R \
        --pondir {OUTPUT}/{wildcards.pon}dir \
        --ref {ref_r} \
        --input {input.vcf} \
        --output {params.vcf} \
        --plotdir {OUTPUT}/plots \
        --scriptdir {SCRIPTS}
        
        rename bgz gz *bgz*
        """



###############################################################################
### OLD STUFF
###############################################################################
rule gridss_normals_preprocessing:
    input:
        # bam=expand(INPUT_BAM+"/{normal}.recalibrated.bam", normal=NORMALS),
        # bam=expand("{normal}", normal=NORMALS),
        # bam=lambda wildcards: expand("{normal}", normal=config[wildcards.pon]),
        bam=NORMAL_BAM+"/{normal}.bam",
        blacklist=res+"/ENCODE_DAC_blacklist/ENCFF419RSJ.bed"
    output:
        bam=OUTPUT+"/{pon}_intermediate_results/{normal}.bam.gridss.working/{normal}.bam.sv.bam"
    threads: 8
    shell:
        """
        mkdir -p {OUTPUT}/{wildcards.pon}_intermediate_results
        ulimit -n 4096

        {SCRIPTS}/gridss.sh \
        -t {threads} \
        -b {input.blacklist} \
        -r {ref} \
        -w {OUTPUT}/{wildcards.pon}_intermediate_results \
        -o {OUTPUT}/{wildcards.pon}.{wildcards.normal}.tmp.vcf \
        -a {OUTPUT}/{wildcards.pon}.gridss.{PROJECT}.assembly.bam \
        -j {SCRIPTS}/gridss-{GRIDSS_VERSION}-gridss-jar-with-dependencies.jar \
        --jvmheap {mem}g \
        --steps PreProcess \
        {input.bam}
        """
        # export _JAVA_OPTIONS="-Dgridss.defensiveGC=true"



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
        {input.tumor} {input.normal}
        """





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