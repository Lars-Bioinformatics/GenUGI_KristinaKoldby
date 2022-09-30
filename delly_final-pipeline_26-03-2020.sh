#!/bin/sh
#
#SBATCH --account sdukoldby_slim      # account
#SBATCH --nodes 1                 # number of nodes
#SBATCH --time 4:00:00            # max time (HH:MM:SS)

####################################################################
############# FINAL PIPELINE USED FOR FURTHER ANALYSES #############
####################################################################

##################################
######### Variant calling ########
##################################

# delly call -g /work/sdukoldby/resources/hg38/Homo_sapiens_assembly38.fasta -q 20 -x delly/human.hg38.excl.tsv -o delly/G45-ECV2-29_raw_delly.bcf G45-ECV2-29-biopsi-C1_truseq-nano-genome_HGL2LDSXX_S21.recalibrated.bam G45-ECV2-29-blod_truseq-nano-genome_HGL2LDSXX_S30.recalibrated.bam
# delly call -g /work/sdukoldby/resources/hg38/Homo_sapiens_assembly38.fasta -q 20 -x delly/human.hg38.excl.tsv -o delly/G45-ECV2-31_raw_delly.bcf G45-ECV2-31-biopsi-F1_truseq-nano-genome_HGL2LDSXX_S22.recalibrated.bam G45-ECV2-31-blod_truseq-nano-genome_HGL2LDSXX_S27.recalibrated.bam
# delly call -g /work/sdukoldby/resources/hg38/Homo_sapiens_assembly38.fasta -q 20 -x delly/human.hg38.excl.tsv -o delly/G45-ECV2-35_raw_delly.bcf G45-ECV2-35-biopsi-G1_truseq-nano-genome_HGL2LDSXX_S23.recalibrated.bam G45-ECV2-35-blod_truseq-nano-genome_HGL2LDSXX_S28.recalibrated.bam
# delly call -g /work/sdukoldby/resources/hg38/Homo_sapiens_assembly38.fasta -q 20 -x delly/human.hg38.excl.tsv -o delly/G45-ECV2-4_raw_delly.bcf G45-ECV2-4-biopsi-H2_truseq-nano-genome_HGL2LDSXX_S24.recalibrated.bam G45-ECV2-4-blod_truseq-nano-genome_HGL2LDSXX_S26.recalibrated.bam
# delly call -g /work/sdukoldby/resources/hg38/Homo_sapiens_assembly38.fasta -q 20 -x delly/human.hg38.excl.tsv -o delly/G45-ECV2-8_raw_delly.bcf G45-ECV2-8-biopsi-I2_truseq-nano-genome_HGL2LDSXX_S25.recalibrated.bam G45-ECV2-8-blod_truseq-nano-genome_HGL2LDSXX_S29.recalibrated.bam

##################################################
########### Somatic-germline filtering ###########
##################################################

# delly filter -f somatic -o delly/G45-ECV2-29_delly_somatic.bcf -s delly/samples.tsv delly/G45-ECV2-29_raw_delly.bcf
# delly filter -f somatic -o delly/G45-ECV2-31_delly_somatic.bcf -s delly/samples.tsv delly/G45-ECV2-31_raw_delly.bcf
# delly filter -f somatic -o delly/G45-ECV2-35_delly_somatic.bcf -s delly/samples.tsv delly/G45-ECV2-35_raw_delly.bcf
# delly filter -f somatic -o delly/G45-ECV2-4_delly_somatic.bcf -s delly/samples.tsv delly/G45-ECV2-4_raw_delly.bcf
# delly filter -f somatic -o delly/G45-ECV2-8_delly_somatic.bcf -s delly/samples.tsv delly/G45-ECV2-8_raw_delly.bcf

############ Bcf to vcf conversion #########

# bcftools view delly/G45-ECV2-29_delly_somatic.bcf > delly/G45-ECV2-29_delly_somatic.vcf
# bcftools view delly/G45-ECV2-31_delly_somatic.bcf > delly/G45-ECV2-31_delly_somatic.vcf
# bcftools view delly/G45-ECV2-35_delly_somatic.bcf > delly/G45-ECV2-35_delly_somatic.vcf
# bcftools view delly/G45-ECV2-4_delly_somatic.bcf > delly/G45-ECV2-4_delly_somatic.vcf
# bcftools view delly/G45-ECV2-8_delly_somatic.bcf > delly/G45-ECV2-8_delly_somatic.vcf

##############################################
##### Genotype tumor SV sites in normals #####
##############################################

# delly call \
# -g /work/sdukoldby/resources/hg38/Homo_sapiens_assembly38.fasta \
# -x delly/human.hg38.excl.tsv \
# -v delly/G45-ECV2-35_delly_somatic.bcf \
# -o delly/G45-ECV2-35-delly_somatic-sites_pon.bcf \
# G45-ECV2-35-biopsi-G1_truseq-nano-genome_HGL2LDSXX_S23.recalibrated.bam \
# G45-ECV2-31-blod_truseq-nano-genome_HGL2LDSXX_S27.recalibrated.bam \
# G45-ECV2-29-blod_truseq-nano-genome_HGL2LDSXX_S30.recalibrated.bam \
# G45-ECV2-35-blod_truseq-nano-genome_HGL2LDSXX_S28.recalibrated.bam \
# G45-ECV2-4-blod_truseq-nano-genome_HGL2LDSXX_S26.recalibrated.bam \
# G45-ECV2-8-blod_truseq-nano-genome_HGL2LDSXX_S29.recalibrated.bam
#
# -v delly/G45-ECV2-4_delly_somatic.bcf \
# -o delly/G45-ECV2-4-delly_somatic-sites_pon.bcf \
# G45-ECV2-4-biopsi-H2_truseq-nano-genome_HGL2LDSXX_S24.recalibrated.bam \

# -v delly/G45-ECV2-8_delly_somatic.bcf \
# -o delly/G45-ECV2-8-delly_somatic-sites_pon.bcf \
# G45-ECV2-8-biopsi-I2_truseq-nano-genome_HGL2LDSXX_S25.recalibrated.bam \

# -v delly/G45-ECV2-29_delly_somatic.bcf \
# -o delly/G45-ECV2-29-delly_somatic-sites_pon.bcf \
# G45-ECV2-29-biopsi-C1_truseq-nano-genome_HGL2LDSXX_S21.recalibrated.bam \

# -v delly/G45-ECV2-31_delly_somatic.bcf \
# -o delly/G45-ECV2-31-delly_somatic-sites_pon.bcf \
# G45-ECV2-31-biopsi-F1_truseq-nano-genome_HGL2LDSXX_S22.recalibrated.bam \


#############################################################
### Post-filter for somatic SVs using all control samples ###
#############################################################

# delly filter -f somatic -o delly/G45-ECV2-4-delly_somatic_pon-filtered.bcf -s delly/samples.tsv delly/G45-ECV2-4-delly_somatic-sites_pon.bcf
delly filter -f somatic -o delly/G45-ECV2-8-delly_somatic_pon-filtered.bcf -s delly/samples.tsv delly/G45-ECV2-8-delly_somatic-sites_pon.bcf
delly filter -f somatic -o delly/G45-ECV2-29-delly_somatic_pon-filtered.bcf -s delly/samples.tsv delly/G45-ECV2-29-delly_somatic-sites_pon.bcf
delly filter -f somatic -o delly/G45-ECV2-31-delly_somatic_pon-filtered.bcf -s delly/samples.tsv delly/G45-ECV2-31-delly_somatic-sites_pon.bcf
delly filter -f somatic -o delly/G45-ECV2-35-delly_somatic_pon-filtered.bcf -s delly/samples.tsv delly/G45-ECV2-35-delly_somatic-sites_pon.bcf

###########################
### Convert bcf to vcf ####
###########################

# bcftools view delly/G45-ECV2-4-delly_somatic_pon-filtered.bcf > delly/G45-ECV2-4-delly_somatic_pon-filtered.vcf
bcftools view delly/G45-ECV2-8-delly_somatic_pon-filtered.bcf > delly/G45-ECV2-8-delly_somatic_pon-filtered.vcf
bcftools view delly/G45-ECV2-29-delly_somatic_pon-filtered.bcf > delly/G45-ECV2-29-delly_somatic_pon-filtered.vcf
bcftools view delly/G45-ECV2-31-delly_somatic_pon-filtered.bcf > delly/G45-ECV2-31-delly_somatic_pon-filtered.vcf
bcftools view delly/G45-ECV2-35-delly_somatic_pon-filtered.bcf > delly/G45-ECV2-35-delly_somatic_pon-filtered.vcf


#####################################################################################################################
######################### TESTED BUT NOT USED FOR FURTHER ANALYSES ##################################################
#####################################################################################################################


#########################
##### Joint calling #####
#########################

# delly call \
# -g /work/sdukoldby/resources/hg38/Homo_sapiens_assembly38.fasta \
# -q 20 \
# -x delly/human.hg38.excl.tsv \
# -o delly/G45-geno_raw_delly.bcf \
# G45-ECV2-31-biopsi-F1_truseq-nano-genome_HGL2LDSXX_S22.recalibrated.bam \
# G45-ECV2-29-biopsi-C1_truseq-nano-genome_HGL2LDSXX_S21.recalibrated.bam \
# G45-ECV2-35-biopsi-G1_truseq-nano-genome_HGL2LDSXX_S23.recalibrated.bam \
# G45-ECV2-4-biopsi-H2_truseq-nano-genome_HGL2LDSXX_S24.recalibrated.bam \
# G45-ECV2-8-biopsi-I2_truseq-nano-genome_HGL2LDSXX_S25.recalibrated.bam \
# G45-ECV2-31-blod_truseq-nano-genome_HGL2LDSXX_S27.recalibrated.bam \
# G45-ECV2-29-blod_truseq-nano-genome_HGL2LDSXX_S30.recalibrated.bam \
# G45-ECV2-35-blod_truseq-nano-genome_HGL2LDSXX_S28.recalibrated.bam \
# G45-ECV2-4-blod_truseq-nano-genome_HGL2LDSXX_S26.recalibrated.bam \
# G45-ECV2-8-blod_truseq-nano-genome_HGL2LDSXX_S29.recalibrated.bam


########### Somatic-germline filtering ###########

# delly filter -f somatic -o delly/G45-ECV2-joint_delly_somatic.bcf -s delly/samples.tsv delly/G45-geno_raw_delly.bcf
# delly filter -p -f somatic -o delly/G45-ECV2-joint_delly_somatic_pass.bcf -s delly/samples.tsv delly/G45-geno_raw_delly.bcf


############ Bcf to vcf conversion #########

# bcftools view delly/G45-ECV2-joint_delly_somatic.bcf > delly/G45-ECV2-joint_delly_somatic.vcf
# bcftools view delly/G45-ECV2-joint_delly_somatic_pass.bcf > delly/G45-ECV2-joint_delly_somatic_pass.vcf


####################################
##### Delly "Panel of Normals" #####
####################################

# delly call \
# -g /work/sdukoldby/resources/hg38/Homo_sapiens_assembly38.fasta \
# -q 20 \
# -x delly/human.hg38.excl.tsv \
# -o delly/G45-normals_raw_delly.bcf \
# G45-ECV2-31-blod_truseq-nano-genome_HGL2LDSXX_S27.recalibrated.bam \
# G45-ECV2-29-blod_truseq-nano-genome_HGL2LDSXX_S30.recalibrated.bam \
# G45-ECV2-35-blod_truseq-nano-genome_HGL2LDSXX_S28.recalibrated.bam \
# G45-ECV2-4-blod_truseq-nano-genome_HGL2LDSXX_S26.recalibrated.bam \
# G45-ECV2-8-blod_truseq-nano-genome_HGL2LDSXX_S29.recalibrated.bam
