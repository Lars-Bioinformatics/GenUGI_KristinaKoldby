__title__ = "Pipeline for Somatic JOINT Variant Calling with Mutect2 - G88-2019 project"
__author__ = "Lars Andersen <larsmew@gmail.com> and Kristina Koldby"
__date__ = "17/03/2020"
__version__ = "1.0"

import time, os, sys, glob

#########################################################
####                       Input                     ####
#########################################################
# Matched tumour-normal samples information
# configfile: "samples.yaml"
# configfile: "somatic_matched_samples_oneSample.yaml"

# Explicit paths for external input files
ref = "/work/sdukoldby/resources/hg38/Homo_sapiens_assembly38.fasta"
gnomad = "/work/sdukoldby/resources/hg38/af-only-gnomad.hg38.vcf.gz"
interval_list = "/work/sdukoldby/resources/hg38/MedExome_hg38_capture_targets.interval_list"
dbsnp = "/work/sdukoldby/resources/hg38/Homo_sapiens_assembly38.dbsnp138.vcf"


SAMPLES = ["G88"]
TUMORS = ["G88-F9B1-1838404-01-05","G88-F9B1-1838404-01-13","G88-F9B1-92-3359-Hud","G88-F9B1-RH07-8852"]
# NORMALS = ["G88-F9B1-Blod"]
NORMALS = ["G88-F9B1-92-3359-Hud"]

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
# mem = "-Xmx12g" # login nodes - be careful not running too many jobs at once!
mem = "-Xmx24g" # slim nodes
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
        expand(output_germline+"{sample}_nimblegen-medexome_germline_haplotypeCaller.vcf", sample=NORMALS),


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
		vcf=output_germline+"{sample}_nimblegen-medexome_germline_haplotypeCaller.vcf"
	shell:
		"""
		gatk --java-options {mem} HaplotypeCaller \
		-R={ref} \
		-I={input.bam} \
		-O={output.vcf} \
		--dbsnp={dbsnp} \
		-L={interval_list} \
		-A Coverage \
		-A DepthPerAlleleBySample \
		-A BaseQualityRankSumTest
		"""
		# Coverage - Total depth of coverage per sample and over all samples
		# DepthPerAlleleBySample - Unfiltered alternative allele depth (AD)
		# BaseQualityRankSumTest - Rank Sum Test of REF versus ALT base quality scores
        # -ERC=GVCF \
