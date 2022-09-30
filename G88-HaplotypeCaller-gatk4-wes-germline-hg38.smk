__title__ = "Pipeline for Somatic JOINT Variant Calling with Mutect2 - G88-2019 project"
__author__ = "Lars Andersen <larsmew@gmail.com> and Kristina Koldby"
__date__ = "17/03/2020"
__version__ = "1.0"

import time, os, sys, glob

#########################################################
####                       Input                     ####
#########################################################
# Matched tumour-normal samples information
configfile: "samples.yaml"
# configfile: "somatic_matched_samples_oneSample.yaml"

# Explicit paths for external input files
ref = "/work/sdukoldby/resources/hg38/Homo_sapiens_assembly38.fasta"
gnomad = "/work/sdukoldby/resources/hg38/af-only-gnomad.hg38.vcf.gz"
interval_list = "/work/sdukoldby/resources/hg38/MedExome_hg38_capture_targets.interval_list"
dbsnp = "/work/sdukoldby/resources/hg38/Homo_sapiens_assembly38.dbsnp138.vcf"


SAMPLES = ["G88-F9B1"]
TUMORS = ["G88-F9B1-1838404-01-05","G88-F9B1-1838404-01-13","G88-F9B1-92-3359-Hud","G88-F9B1-RH07-8852"]
NORMALS = ["G88-F9B1-Blod"]

print(SAMPLES)
print(TUMORS)
print(NORMALS)
# sys.exit()

#########################################################
####                      Output                     ####
#########################################################
log_file = "log_file_germline.txt"
output_germline = "germline_variants/"


#########################################################
####                       Setup                     ####
#########################################################
# Timing
totim = time.time()
timeFormat = "%Y_%m_%d:%X" # year, month, day, time H:M:S

# Memory
mem = "-Xmx12g" # login nodes - be careful not running too many jobs at once!
# mem = "-Xmx24g" # slim nodes
# mem = "-Xmx32g" # Fat nodes
# mem = "-Xmx64g"

onstart:
    shell("mkdir -p "+output_germline)


#########################################################
####                  Run All Rules                  ####
#########################################################
'''
Rule all
'''

rule all:
    input:
        expand(output_germline+"{fam_name}_germline_variants_haplotypeCaller.vcf.gz", fam_name=SAMPLES),


###############################################################################
### Call germline variants                                                  ###
###############################################################################
'''
HaplotypeCaller
'''
rule HaplotypeCaller:
    input:
        bam="{sample}_nimblegen-medexome_HGVM5DSXX.connor.recalibrated.bam"
    output:
        gvcf="gvcf_files/{sample}_variants.g.vcf.gz"
    shell:
        """
        gatk --java-options {mem} HaplotypeCaller \
        -R {ref} \
        -I {input.bam} \
        -O {output.gvcf} \
        -ERC GVCF \
        --dbsnp {dbsnp} \
        -L {interval_list} \
        -A Coverage \
        -A DepthPerAlleleBySample \
        -A BaseQualityRankSumTest
        """
        # Coverage - Total depth of coverage per sample and over all samples
        # DepthPerAlleleBySample - Unfiltered alternative allele depth (AD)
        # BaseQualityRankSumTest - Rank Sum Test of REF versus ALT base quality scores

'''
Create single multi-sample g.vcf file    
'''
rule CombineGVCFs:
    input:
        gvcfs=lambda wildcards: expand("gvcf_files/{sample}_variants.g.vcf.gz", sample=config[wildcards.fam_name]["all"])
    output:
        gvcf="gvcf_files/{fam_name}_variants_combined.g.vcf.gz"
    params:
        gvcfs=lambda wildcards: expand("-V gvcf_files/{sample}_variants.g.vcf.gz", sample=config[wildcards.fam_name]["all"]),
    shell:
        """
        gatk --java-options {mem} CombineGVCFs \
        -R {ref} \
        {params.gvcfs} \
        -O {output.gvcf}
        """

'''
GenotypeGVCFs
'''
rule GenotypeGVCFs:
    input:
        gvcf="gvcf_files/{fam_name}_variants_combined.g.vcf.gz"
    output:
        vcf=output_germline+"{fam_name}_germline_variants_haplotypeCaller.vcf.gz"
    shell:
        """
        gatk --java-options {mem} GenotypeGVCFs \
        -R {ref} \
        -V {input.gvcf} \
        --dbsnp {dbsnp} \
        -L {interval_list} \
        -O {output} \
        """
