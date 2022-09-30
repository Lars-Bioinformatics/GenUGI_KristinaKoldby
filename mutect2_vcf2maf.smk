__title__ = "vcf2maf"
__author__ = "Lars Andersen <larsmew@gmail.com>"
__date__ = "02/06/2021"
__version__ = "1.0"

import time

#########################################################
####                       Input                     ####
#########################################################
# Matched tumour-normal samples information
configfile: "/work/Data/samples-pairwise.yaml"
# configfile: "somatic_matched_samples_oneSample.yaml"

# Explicit paths for external input files
ref = "/work/sduvarcall/resources/hg38/Homo_sapiens_assembly38.fasta"
vep_data = "/work/sduvarcall/ensembl_vep/vep"
vep_path = "/work/miniconda3/envs/vep/bin"


# ref_build = "GRCh37"
ref_build = "GRCh38"


#########################################################
####                      Output                     ####
#########################################################
output_mutect2 = "/work/Data/Connor/mutect2_vcf/"
output_varscan2 = "/work/Data/Connor/varscan_somatic_joint_calling/"

#########################################################
####               Sample Information                ####
#########################################################
# Mutect2 sample information
PAIR = [pair for pair in config]

# Varscan2 sample information
VARSCAN_SAMPLES = glob_wildcards(output_varscan2+"{sample}_cns_varscan2.vcf.gz")
print(VARSCAN_SAMPLES)

#########################################################
####                       Setup                     ####
#########################################################
# Timing
totim = time.time()
timeFormat = "%Y_%m_%d:%X" # year, month, day, time H:M:S

# Memory
mem = "-Xmx12g"

#########################################################
####  Define workflow start, stop and error actions  ####
#########################################################
onstart:
    # shell("mkdir -p "+output_mutect2)
    shell("mkdir -p "+output_mutect2+"mutect2_maf/")
    shell("mkdir -p "+output_mutect2+"mutect2_maf_PASSonly/")
    shell("mkdir -p "+output_varscan2+"varscan2_maf/")
    shell("mkdir -p "+output_varscan2+"varscan2_maf_somatic_pon-filtered/")


#########################################################
####                  Run All Rules                  ####
#########################################################
'''
Rule all
'''
rule all_pairs:
    input:
        [expand(output_mutect2+"mutect2_maf/{tumor}_vs_{normal}_somatic_mutect2_filterFlag.maf",
            normal=config[fam]["normal"],
            tumor=config[fam]["tumor"]) for fam in PAIR],
        [expand(output_mutect2+"mutect2_maf_PASSonly/{tumor}_vs_{normal}_somatic_mutect2_filterFlag_norm_PASSonly.maf",
            normal=config[fam]["normal"],
            tumor=config[fam]["tumor"]) for fam in PAIR],


#########################################################
####                  Run All Rules                  ####
#########################################################
rule vcf2maf_mutect2:
    input:
        vcf=output_mutect2+"{tumor}_vs_{normal}_somatic_mutect2_filterFlag.vcf.gz"
    output:
        maf=output_mutect2+"mutect2_maf/{tumor}_vs_{normal}_somatic_mutect2_filterFlag_norm.maf",
        vep_vcf=output_mutect2+"{tumor}_vs_{normal}_somatic_mutect2_filterFlag_norm.vep.vcf",
        vcf=temp(output_mutect2+"{tumor}_vs_{normal}_somatic_mutect2_filterFlag_norm.vcf")
    #params:
    #    vcf=temp(output_mutect2+"{tumor}_vs_{normal}_somatic_mutect2_filterFlag.vcf")
    shell:
        """
        bgzip -d -c {input.vcf} > {output.vcf}
        vcf2maf.pl \
        --input-vcf {output.vcf} \
        --output-maf {output.maf} \
        --vep-data {vep_data} \
        --vep-path {vep_path} \
        --ncbi-build {ref_build} \
        --ref-fasta {ref} \
        --tumor-id {wildcards.tumor} \
        --normal-id {wildcards.normal} \
        --vcf-tumor-id {wildcards.tumor} \
        --vcf-normal-id {wildcards.normal}
        """

rule vcf2maf_mutect2_PASS:
    input:
        vcf=output_mutect2+"{tumor}_vs_{normal}_somatic_mutect2_filterFlag_PASSonly.vcf.gz"
    output:
        maf=output_mutect2+"mutect2_maf_PASSonly/{tumor}_vs_{normal}_somatic_mutect2_filterFlag_norm_PASSonly.maf",
        vep_vcf=output_mutect2+"{tumor}_vs_{normal}_somatic_mutect2_filterFlag_norm_PASSonly.vep.vcf",
        vcf=temp(output_mutect2+"{tumor}_vs_{normal}_somatic_mutect2_filterFlag_norm_PASSonly.vcf")
    shell:
        """
        bgzip -d -c {input.vcf} > {output.vcf}
        vcf2maf.pl \
        --input-vcf {output.vcf} \
        --output-maf {output.maf} \
        --vep-data {vep_data} \
        --vep-path {vep_path} \
        --ncbi-build {ref_build} \
        --ref-fasta {ref} \
        --tumor-id {wildcards.tumor} \
        --normal-id {wildcards.normal} \
        --vcf-tumor-id {wildcards.tumor} \
        --vcf-normal-id {wildcards.normal}
        """
