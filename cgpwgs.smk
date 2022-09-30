__author__ = "Lars Andersen <larsmew@gmail.com> & Kristina Koldby"
__date__ = "24/01/2020"
__version__ = "2.0"

import time, os, sys

#SAMPLES, = glob_wildcards("{sample}_R1.fastq.gz")
totim = time.time()
timeFormat = time.strftime("%Y%m%d-%H%M%S") # year, month, day, time H:M:S

WORK = "/work/sdukoldby/G45-2016-genugi/wgs_bam/"
RESOURCES = "/work/sdukoldby/resources/hg38"

normal_sample = "ECV2-29_blod"
tumor_sample = "ECV2-29_biopsi-C1"

# print(normal_sample)
# print(tumor_sample)

# ploidy = 3.2
# purity = 0.35

'''
Create log file, containing:
	- Programs used and their version info
	- Start and stop time of complete workflow
	- sample(s)
'''
onstart:
	shell("echo 'Started execution of pipeline:' $(date +'%Y-%m-%d %H:%M:%S')")

onsuccess:
	fiTime = 'Total time:', str((time.time()-totim) / 60), 'minutes'
	shell("echo 'Finished execution on' $(date +'%Y-%m-%d %H:%M:%S')")
	shell("echo {fiTime}")

onerror:
	fiTime = 'Total time:', str((time.time()-totim) / 60), 'minutes'
	shell("echo 'Finished execution on' $(date +'%Y-%m-%d %H:%M:%S')")
	shell("echo {fiTime}")
	shell("echo 'ERROR OCCURED, PLEASE REFER TO SLURM LOGFILE FOR DETAILS'")


'''
Rule all
'''
rule all:
	input:
		# BLOCK 6
		expand(WORK+"/wgs_out/{tumor}_vs_{normal}/timings/WGS_{tumor}_vs_{normal}.time.verify_MT", normal=normal_sample, tumor=tumor_sample),
		expand(WORK+"/wgs_out/{tumor}_vs_{normal}/timings/WGS_{tumor}_vs_{normal}.time.CaVEMan_annot", normal=normal_sample, tumor=tumor_sample),
		expand(WORK+"/wgs_out/{tumor}_vs_{normal}/timings/WGS_{tumor}_vs_{normal}.time.CaVEMan_annot_germline", normal=normal_sample, tumor=tumor_sample),
		expand(WORK+"/wgs_out/{tumor}_vs_{normal}/timings/WGS_{tumor}_vs_{normal}.time.BRASS", normal=normal_sample, tumor=tumor_sample),
		# BLOCK 3
		# expand(WORK+"/wgs_out/{tumor}_vs_{normal}/timings/WGS_{tumor}_vs_{normal}.time.ascat", normal=normal_sample, tumor=tumor_sample),
		# expand(WORK+"/wgs_out/{tumor}_vs_{normal}/timings/WGS_{tumor}_vs_{normal}.time.BRASS_cover", normal=normal_sample, tumor=tumor_sample),
		# expand(WORK+"/wgs_out/{tumor}_vs_{normal}/timings/WGS_{tumor}_vs_{normal}.time.BRASS_input", normal=normal_sample, tumor=tumor_sample),
		# Manual execution ascat with purity and ploidy
		#expand(WORK+"/wgs_out/{tumor}_vs_{normal}/timings/WGS_{tumor}_vs_{normal}.time.ascat_pu_pi", normal=normal_sample, tumor=tumor_sample),
		#expand(WORK+"/wgs_out/{tumor}_vs_{normal}/timings/WGS_{tumor}_vs_{normal}.time.ascat", normal=normal_sample, tumor=tumor_sample),
		#expand(WORK+"/wgs_out/{tumor}_vs_{normal}/timings/WGS_{tumor}_vs_{normal}.time.ascat_test1", normal=normal_sample, tumor=tumor_sample),
		# Copy number - ascat and battenberg
		# expand(WORK+"/wgs_out/{tumor}_vs_{normal}/timings/WGS_{tumor}_vs_{normal}.time.ascat", normal=normal_sample, tumor=tumor_sample),
# 		expand(WORK+"/wgs_out/{tumor}_vs_{normal}/timings/WGS_{tumor}_vs_{normal}.time.battenberg", normal=normal_sample, tumor=tumor_sample)




###############################################################################
#### Parallel Block 1
#### CaVEMan setup
#### BB splitlocifiles
#### Genotype Check
#### VerifyBam Normal
###############################################################################

rule caveman_setup:
	input:
		normal=WORK+"/bam/{normal}.bam",
		tumor=WORK+"/bam/{tumor}.bam"
	output:
		res=WORK+"/wgs_out/{tumor}_vs_{normal}/WGS_{tumor}_vs_{normal}/genotyped/result.json"
	shell:
		"""

		mkdir -p {WORK}/wgs_out/{wildcards.tumor}_vs_{wildcards.normal}; \
		mkdir -p {WORK}/wgs_out/{wildcards.tumor}_vs_{wildcards.normal}/timings; \
		udocker run \
		cgpwgs \
		bash -c '/usr/bin/time -v compareBamGenotypes.pl  \
		-o /work/wgs_out/{wildcards.tumor}_vs_{wildcards.normal}/WGS_{wildcards.tumor}_vs_{wildcards.normal}/genotyped  \
		-j /work/wgs_out/{wildcards.tumor}_vs_{wildcards.normal}/WGS_{wildcards.tumor}_vs_{wildcards.normal}/genotyped/result.json \
		-nb /work/bam_out/{wildcards.normal}/{wildcards.normal}.bam \
		-tb /work/bam_out/{wildcards.tumor}/{wildcards.tumor}.bam \
		-s /work/reference_files/general.tsv  \
		-g /work/reference_files/gender.tsv \
		>& /work/wgs_out/{wildcards.tumor}_vs_{wildcards.normal}/timings/WGS_{wildcards.tumor}_vs_{wildcards.normal}.time.geno ;
		echo \''WRAPPER_EXIT: '\' $?'
		"""


rule bb_split_loci_files:
	input:
		normal=WORK+"/bam/{normal}.bam",
		tumor=WORK+"/bam/{tumor}.bam"
	output:
		#lst=WORK+"/wgs_out/{tumor}_vs_{normal}/WGS_{tumor}_vs_{normal}/caveman/splitList"
		timing=WORK+"/wgs_out/{tumor}_vs_{normal}/timings/WGS_{tumor}_vs_{normal}.time.CaVEMan_setup"
	threads: 24
	shell:
		"""
		mkdir -p {WORK}/wgs_out/{wildcards.tumor}_vs_{wildcards.normal}; \
		mkdir -p {WORK}/wgs_out/{wildcards.tumor}_vs_{wildcards.normal}/timings; \
		mkdir -p {WORK}/wgs_out/{wildcards.tumor}_vs_{wildcards.normal}/tmp; \
		cp {WORK}/reference_files/caveman/bed_dummies_for_setup/*bed {WORK}/wgs_out/{wildcards.tumor}_vs_{wildcards.normal}/tmp/; \
		cp {WORK}/reference_files/caveman/flag.vcf.config.WGS.ini {WORK}/wgs_out/{wildcards.tumor}_vs_{wildcards.normal}/; \
		python2 ~/udocker run \
		--workdir=/  \
		--user=laran  \
		--volume={WORK}:/work  \
		cgpwgs \
		bash -c '/usr/bin/time -v caveman.pl \
		-r /work/reference_files/genome.fa.fai \
		-ig /work/reference_files/caveman/HiDepth.tsv\
		-b /work/reference_files/caveman/flagging \
		-f /work/reference_files/caveman/flagging/flag.to.vcf.convert.ini \
		-c /work/wgs_out/{wildcards.tumor}_vs_{wildcards.normal}/flag.vcf.config.WGS.ini \
		-ab /work/reference_files/vagrent \
		-u /work/reference_files/caveman \
		-s '\''Human'\'' \
		-sa NCBI37 \
		-t {threads} \
		-st WGS \
		-td 5 \
		-nd 2 \
		-nb /work/bam_out/{wildcards.normal}/{wildcards.normal}.bam \
		-tb /work/bam_out/{wildcards.tumor}/{wildcards.tumor}.bam \
		-tc /work/wgs_out/{wildcards.tumor}_vs_{wildcards.normal}/tmp/tum.cn.bed \
		-nc /work/wgs_out/{wildcards.tumor}_vs_{wildcards.normal}/tmp/norm.cn.bed \
		-e 350000 \
		-o /work/wgs_out/{wildcards.tumor}_vs_{wildcards.normal}/WGS_{wildcards.tumor}_vs_{wildcards.normal}/caveman \
		-x NC_007605,hs37d5,GL% \
		-p setup \
		>& /work/wgs_out/{wildcards.tumor}_vs_{wildcards.normal}/timings/WGS_{wildcards.tumor}_vs_{wildcards.normal}.time.CaVEMan_setup ;
		echo \''WRAPPER_EXIT: '\' $?'
		"""

rule genotype_check:
	input:
		normal=WORK+"/bam_out/{normal}/{normal}.bam",
		tumor=WORK+"/bam_out/{tumor}/{tumor}.bam"
	output:
		timing=WORK+"/wgs_out/{tumor}_vs_{normal}/timings/WGS_{tumor}_vs_{normal}.time.splitlocifiles"
	threads: 24
	shell:
		"""
		mkdir -p {WORK}/wgs_out/{wildcards.tumor}_vs_{wildcards.normal}; \
		mkdir -p {WORK}/wgs_out/{wildcards.tumor}_vs_{wildcards.normal}/timings; \
		python2 ~/udocker run \
		--workdir=/  \
		--user=laran  \
		--volume={WORK}:/work  \
		cgpwgs \
		bash -c '/usr/bin/time -v battenberg.pl \
		-o /work/wgs_out/{wildcards.tumor}_vs_{wildcards.normal}/WGS_{wildcards.tumor}_vs_{wildcards.normal}/battenberg \
		-u /work/reference_files/battenberg/1000genomesloci \
		-e /work/reference_files/battenberg/impute/impute_info.txt \
		-c /work/reference_files/battenberg/probloci.txt \
		-r /work/reference_files/genome.fa.fai \
		-ig /work/reference_files/battenberg/ignore_contigs.txt \
		-ge XX \
		-tb /work/bam_out/{wildcards.tumor}/{wildcards.tumor}.bam \
		-nb /work/bam_out/{wildcards.normal}/{wildcards.normal}.bam \
		-p splitlocifiles \
		-nl 50 \
		-t {threads} \
		>& /work/wgs_out/{wildcards.tumor}_vs_{wildcards.normal}/timings/WGS_{wildcards.tumor}_vs_{wildcards.normal}.time.splitlocifiles ;
		echo \''WRAPPER_EXIT: '\' $?'
		"""

rule verifybam_normal:
	input:
		normal=WORK+"/bam_out/{normal}/{normal}.bam",
		tumor=WORK+"/bam_out/{tumor}/{tumor}.bam"
	output:
		timing=WORK+"/wgs_out/{tumor}_vs_{normal}/timings/WGS_{tumor}_vs_{normal}.time.verify_WT"
	threads: 24
	shell:
		"""
		mkdir -p {WORK}/wgs_out/{wildcards.tumor}_vs_{wildcards.normal}; \
		mkdir -p {WORK}/wgs_out/{wildcards.tumor}_vs_{wildcards.normal}/timings; \
		python2 ~/udocker run \
		--workdir=/  \
		--user=laran  \
		--volume={WORK}:/work  \
		cgpwgs \
		bash -c '/usr/bin/time -v verifyBamHomChk.pl \
		-d 25 \
		-o /work/wgs_out/{wildcards.tumor}_vs_{wildcards.normal}/WGS_{wildcards.normal}/contamination \
		-b /work/bam_out/{wildcards.normal}/{wildcards.normal}.bam \
		-t {threads} \
		-j /work/wgs_out/{wildcards.tumor}_vs_{wildcards.normal}/WGS_{wildcards.normal}/contamination/result.json \
		-s /work/reference_files/verifyBamID_snps.vcf.gz \
		>& /work/wgs_out/{wildcards.tumor}_vs_{wildcards.normal}/timings/WGS_{wildcards.tumor}_vs_{wildcards.normal}.time.verify_WT ;
		echo \''WRAPPER_EXIT: '\' $?'
		"""


###############################################################################
#### Parallel Block 2
#### BB alleleCount
#### CaVEMan split
###############################################################################
rule caveman_split:
	input:
		normal=WORK+"/bam_out/{normal}/{normal}.bam",
		tumor=WORK+"/bam_out/{tumor}/{tumor}.bam",
		res1_1=WORK+"/wgs_out/{tumor}_vs_{normal}/WGS_{tumor}_vs_{normal}/genotyped/result.json",
		res1_2=WORK+"/wgs_out/{tumor}_vs_{normal}/timings/WGS_{tumor}_vs_{normal}.time.CaVEMan_setup",
		res1_3=WORK+"/wgs_out/{tumor}_vs_{normal}/timings/WGS_{tumor}_vs_{normal}.time.splitlocifiles",
		res1_4=WORK+"/wgs_out/{tumor}_vs_{normal}/timings/WGS_{tumor}_vs_{normal}.time.verify_WT"
	output:
		timing=WORK+"/wgs_out/{tumor}_vs_{normal}/timings/WGS_{tumor}_vs_{normal}.time.CaVEMan_split"
	threads: 24
	shell:
		"""
		python2 ~/udocker run \
		--workdir=/  \
		--user=laran  \
		--volume={WORK}:/work  \
		cgpwgs \
		bash -c '/usr/bin/time -v caveman.pl \
		-r /work/reference_files/genome.fa.fai \
		-ig /work/reference_files/caveman/HiDepth.tsv \
		-b /work/reference_files/caveman/flagging \
		-ab /work/reference_files/vagrent \
		-u /work/reference_files/caveman \
		-s '\''Human'\'' \
		-sa NCBI37 \
		-t {threads} \
		-st WGS \
		-tc /work/wgs_out/{wildcards.tumor}_vs_{wildcards.normal}/tmp/tum.cn.bed \
		-nc /work/wgs_out/{wildcards.tumor}_vs_{wildcards.normal}/tmp/norm.cn.bed \
		-td 5 \
		-nd 2 \
		-tb /work/bam_out/{wildcards.tumor}/{wildcards.tumor}.bam \
		-nb /work/bam_out/{wildcards.normal}/{wildcards.normal}.bam \
		-c /work/wgs_out/{wildcards.tumor}_vs_{wildcards.normal}/flag.vcf.config.WGS.ini \
		-f /work/reference_files/caveman/flagging/flag.to.vcf.convert.ini \
		-e 350000 \
		-o /work/wgs_out/{wildcards.tumor}_vs_{wildcards.normal}/WGS_{wildcards.tumor}_vs_{wildcards.normal}/caveman \
		-x NC_007605,hs37d5,GL% \
		-p split \
		>& /work/wgs_out/{wildcards.tumor}_vs_{wildcards.normal}/timings/WGS_{wildcards.tumor}_vs_{wildcards.normal}.time.CaVEMan_split ;
		echo \''WRAPPER_EXIT: '\' $?'
		"""

rule bb_allelecount:
	input:
		normal=WORK+"/bam_out/{normal}/{normal}.bam",
		tumor=WORK+"/bam_out/{tumor}/{tumor}.bam",
		res1_1=WORK+"/wgs_out/{tumor}_vs_{normal}/WGS_{tumor}_vs_{normal}/genotyped/result.json",
		res1_2=WORK+"/wgs_out/{tumor}_vs_{normal}/timings/WGS_{tumor}_vs_{normal}.time.CaVEMan_setup",
		res1_3=WORK+"/wgs_out/{tumor}_vs_{normal}/timings/WGS_{tumor}_vs_{normal}.time.splitlocifiles",
		res1_4=WORK+"/wgs_out/{tumor}_vs_{normal}/timings/WGS_{tumor}_vs_{normal}.time.verify_WT"
	output:
		timing=WORK+"/wgs_out/{tumor}_vs_{normal}/timings/WGS_{tumor}_vs_{normal}.time.alleleCount"
	threads: 24
	shell:
		"""
		python2 ~/udocker run \
		--workdir=/  \
		--user=laran  \
		--volume={WORK}:/work  \
		cgpwgs \
		bash -c '/usr/bin/time -v battenberg.pl \
		-o /work/wgs_out/{wildcards.tumor}_vs_{wildcards.normal}/WGS_{wildcards.tumor}_vs_{wildcards.normal}/battenberg \
		-u /work/reference_files/battenberg/1000genomesloci \
		-e /work/reference_files/battenberg/impute/impute_info.txt \
		-c /work/reference_files/battenberg/probloci.txt \
		-r /work/reference_files/genome.fa.fai \
		-ig /work/reference_files/battenberg/ignore_contigs.txt \
		-ge XX \
		-tb /work/bam_out/{wildcards.tumor}/{wildcards.tumor}.bam \
		-nb /work/bam_out/{wildcards.normal}/{wildcards.normal}.bam \
		-p allelecount \
		-nl 50 \
		-t {threads} \
		>& /work/wgs_out/{wildcards.tumor}_vs_{wildcards.normal}/timings/WGS_{wildcards.tumor}_vs_{wildcards.normal}.time.alleleCount ;
		echo \''WRAPPER_EXIT: '\' $?'
		"""

###############################################################################
#### Parallel Block 3
#### ASCAT
#### BRASS_input
#### BRASS_cover
###############################################################################
rule ascat:
	input:
		normal=WORK+"/bam_out/{normal}/{normal}.bam",
		tumor=WORK+"/bam_out/{tumor}/{tumor}.bam",
		res2_1=WORK+"/wgs_out/{tumor}_vs_{normal}/timings/WGS_{tumor}_vs_{normal}.time.CaVEMan_split",
		res2_2=WORK+"/wgs_out/{tumor}_vs_{normal}/timings/WGS_{tumor}_vs_{normal}.time.alleleCount"
	output:
		timing=WORK+"/wgs_out/{tumor}_vs_{normal}/timings/WGS_{tumor}_vs_{normal}.time.ascat"
	threads: 24
	shell:
		"""
		python2 ~/udocker run \
		--workdir=/  \
		--user=laran  \
		--volume={WORK}:/work  \
		cgpwgs \
		bash -c '/usr/bin/time -v ascat.pl \
		-o /work/wgs_out/{wildcards.tumor}_vs_{wildcards.normal}/WGS_{wildcards.tumor}_vs_{wildcards.normal}/ascat \
		-t /work/bam_out/{wildcards.tumor}/{wildcards.tumor}.bam \
		-n /work/bam_out/{wildcards.normal}/{wildcards.normal}.bam \
		-sg /work/reference_files/ascat/SnpGcCorrections.tsv \
		-r /work/reference_files/genome.fa \
		-q 20 \
		-g L \
		-rs '\''Human'\'' \
		-ra NCBI37 \
		-pr WGS \
		-pl ILLUMINA \
		-c {threads} \
		>& /work/wgs_out/{wildcards.tumor}_vs_{wildcards.normal}/timings/WGS_{wildcards.tumor}_vs_{wildcards.normal}.time.ascat ;
		echo \''WRAPPER_EXIT: '\' $?; \
		cut -f 2,3,4,5 -d "," \
		/work/wgs_out/{wildcards.tumor}_vs_{wildcards.normal}/WGS_{wildcards.tumor}_vs_{wildcards.normal}/ascat/{wildcards.tumor}.copynumber.caveman.csv \
		| tr -s "," "\t" > /work/wgs_out/{wildcards.tumor}_vs_{wildcards.normal}/tmp/norm_ascat.cn.bed;
		cut -f 2,3,4,7 -d "," \
		/work/wgs_out/{wildcards.tumor}_vs_{wildcards.normal}/WGS_{wildcards.tumor}_vs_{wildcards.normal}/ascat/{wildcards.tumor}.copynumber.caveman.csv \
		| tr -s "," "\t" > /work/wgs_out/{wildcards.tumor}_vs_{wildcards.normal}/tmp/tum_ascat.cn.bed'
		"""

# Ascat with fixed ploidy and purity - needs to be run manually
rule ascat_pu_pi:
	input:
		normal=WORK+"/bam_out/{normal}/{normal}.bam",
		tumor=WORK+"/bam_out/{tumor}/{tumor}.bam",
		res2_1=WORK+"/wgs_out/{tumor}_vs_{normal}/timings/WGS_{tumor}_vs_{normal}.time.CaVEMan_split",
		res2_2=WORK+"/wgs_out/{tumor}_vs_{normal}/timings/WGS_{tumor}_vs_{normal}.time.alleleCount"
	output:
		timing=WORK+"/wgs_out/{tumor}_vs_{normal}/timings/WGS_{tumor}_vs_{normal}.time.ascat_pu_pi"
	threads: 24
	params:
		purity=purity,
		ploidy=ploidy
	shell:
		"""
		python2 ~/udocker run \
		--workdir=/  \
		--user=laran  \
		--volume={WORK}:/work  \
		cgpwgs \
		bash -c '/usr/bin/time -v ascat.pl \
		-o /work/wgs_out/{wildcards.tumor}_vs_{wildcards.normal}/WGS_{wildcards.tumor}_vs_{wildcards.normal}/ascat_pu_pi \
		-t /work/bam_out/{wildcards.tumor}/{wildcards.tumor}.bam \
		-n /work/bam_out/{wildcards.normal}/{wildcards.normal}.bam \
		-sg /work/reference_files/ascat/SnpGcCorrections.tsv \
		-r /work/reference_files/genome.fa \
		-q 20 \
		-g L \
		-rs '\''Human'\'' \
		-ra NCBI37 \
		-pr WGS \
		-pl ILLUMINA \
		-c {threads} \
		-pi {params.ploidy} \
		-pu {params.purity} \
		>& /work/wgs_out/{wildcards.tumor}_vs_{wildcards.normal}/timings/WGS_{wildcards.tumor}_vs_{wildcards.normal}.time.ascat_pu_pi ;
		echo \''WRAPPER_EXIT: '\' $?; \
		cut -f 2,3,4,5 -d "," \
		/work/wgs_out/{wildcards.tumor}_vs_{wildcards.normal}/WGS_{wildcards.tumor}_vs_{wildcards.normal}/ascat_pu_pi/{wildcards.tumor}.copynumber.caveman.csv \
		| tr -s "," "\t" > /work/wgs_out/{wildcards.tumor}_vs_{wildcards.normal}/tmp/norm_pu_pi.cn.bed;
		cut -f 2,3,4,7 -d "," \
		/work/wgs_out/{wildcards.tumor}_vs_{wildcards.normal}/WGS_{wildcards.tumor}_vs_{wildcards.normal}/ascat_pu_pi/{wildcards.tumor}.copynumber.caveman.csv \
		| tr -s "," "\t" > /work/wgs_out/{wildcards.tumor}_vs_{wildcards.normal}/tmp/tum_pu_pi.cn.bed'
		"""


rule Brass_cover:
	input:
		normal=WORK+"/bam_out/{normal}/{normal}.bam",
		tumor=WORK+"/bam_out/{tumor}/{tumor}.bam",
		res2_1=WORK+"/wgs_out/{tumor}_vs_{normal}/timings/WGS_{tumor}_vs_{normal}.time.CaVEMan_split",
		res2_2=WORK+"/wgs_out/{tumor}_vs_{normal}/timings/WGS_{tumor}_vs_{normal}.time.alleleCount"
	output:
		timing=WORK+"/wgs_out/{tumor}_vs_{normal}/timings/WGS_{tumor}_vs_{normal}.time.BRASS_cover"
	threads: 24
	shell:
		"""
		python2 ~/udocker run \
		--workdir=/  \
		--user=laran  \
		--volume={WORK}:/work  \
		cgpwgs \
		bash -c '/usr/bin/time -v nice -n 10 brass.pl \
		-j 4 -k 4 -c {threads} \
		-d /work/reference_files/brass/HiDepth.bed.gz \
		-f /work/reference_files/brass/brass_np.groups.gz \
		-g /work/reference_files/genome.fa \
		-s '\''Human'\'' -as NCBI37 -pr WGS -pl ILLUMINA \
		-g_cache /work/reference_files/vagrent/vagrent.cache.gz \
		-vi /work/reference_files/brass/viral.genomic.fa.2bit \
		-mi /work/reference_files/brass/all_ncbi_bacteria \
		-b /work/reference_files/brass/500bp_windows.gc.bed.gz \
		-ct /work/reference_files/brass/CentTelo.tsv \
		-cb /work/reference_files/brass/cytoband.txt \
		-t /work/bam_out/{wildcards.tumor}/{wildcards.tumor}.bam \
		-n /work/bam_out/{wildcards.normal}/{wildcards.normal}.bam \
		-o /work/wgs_out/{wildcards.tumor}_vs_{wildcards.normal}/WGS_{wildcards.tumor}_vs_{wildcards.normal}/brass \
		-p cover \
		>& /work/wgs_out/{wildcards.tumor}_vs_{wildcards.normal}/timings/WGS_{wildcards.tumor}_vs_{wildcards.normal}.time.BRASS_cover ;
		echo \''WRAPPER_EXIT: '\' $?'
		"""

rule Brass_input:
	input:
		normal=WORK+"/bam_out/{normal}/{normal}.bam",
		tumor=WORK+"/bam_out/{tumor}/{tumor}.bam",
		res2_1=WORK+"/wgs_out/{tumor}_vs_{normal}/timings/WGS_{tumor}_vs_{normal}.time.CaVEMan_split",
		res2_2=WORK+"/wgs_out/{tumor}_vs_{normal}/timings/WGS_{tumor}_vs_{normal}.time.alleleCount"
	output:
		timing=WORK+"/wgs_out/{tumor}_vs_{normal}/timings/WGS_{tumor}_vs_{normal}.time.BRASS_input"
	threads: 24
	shell:
		"""
		python2 ~/udocker run \
		--workdir=/  \
		--user=laran  \
		--volume={WORK}:/work  \
		cgpwgs \
		bash -c '/usr/bin/time -v brass.pl \
		-j 4 -k 4 -c {threads} \
		-d /work/reference_files/brass/HiDepth.bed.gz \
		-f /work/reference_files/brass/brass_np.groups.gz \
		-g /work/reference_files/genome.fa \
		-s '\''Human'\'' -as NCBI37 -pr WGS -pl ILLUMINA \
		-g_cache /work/reference_files/vagrent/vagrent.cache.gz \
		-vi /work/reference_files/brass/viral.genomic.fa.2bit \
		-mi /work/reference_files/brass/all_ncbi_bacteria \
		-b /work/reference_files/brass/500bp_windows.gc.bed.gz \
		-ct /work/reference_files/brass/CentTelo.tsv \
		-cb /work/reference_files/brass/cytoband.txt \
		-t /work/bam_out/{wildcards.tumor}/{wildcards.tumor}.bam \
		-n /work/bam_out/{wildcards.normal}/{wildcards.normal}.bam \
		-o /work/wgs_out/{wildcards.tumor}_vs_{wildcards.normal}/WGS_{wildcards.tumor}_vs_{wildcards.normal}/brass \
		-p input \
		>& /work/wgs_out/{wildcards.tumor}_vs_{wildcards.normal}/timings/WGS_{wildcards.tumor}_vs_{wildcards.normal}.time.BRASS_input ;
		echo  \''WRAPPER_EXIT: '\' $?'
		"""

rule convert_battenberg_csv:
	input:
		res2_1=WORK+"/wgs_out/{tumor}_vs_{normal}/timings/WGS_{tumor}_vs_{normal}.time.CaVEMan_split",
		res2_2=WORK+"/wgs_out/{tumor}_vs_{normal}/timings/WGS_{tumor}_vs_{normal}.time.alleleCount",
		battenberg=WORK+"/wgs_out/{tumor}_vs_{normal}/timings/WGS_{tumor}_vs_{normal}.time.battenberg"
	output:
		norm_bed=WORK+"/wgs_out/{tumor}_vs_{normal}/tmp/norm_bat.cn.bed",
		tum_bed=WORK+"/wgs_out/{tumor}_vs_{normal}/tmp/tum_bat.cn.bed"
	threads: 24
	shell:
		"""
		cut -f 2,3,4,5 -d "," \
		{WORK}/wgs_out/{wildcards.tumor}_vs_{wildcards.normal}/battenberg/{wildcards.tumor}_battenberg_copynumber.csv \
		| tr -s "," "\t" > {WORK}/wgs_out/{wildcards.tumor}_vs_{wildcards.normal}/tmp/norm_bat.cn.bed; \
		cut -f 2,3,4,7 -d "," \
		{WORK}/wgs_out/{wildcards.tumor}_vs_{wildcards.normal}/battenberg/{wildcards.tumor}_battenberg_copynumber.csv \
		| tr -s "," "\t" > {WORK}/wgs_out/{wildcards.tumor}_vs_{wildcards.normal}/tmp/tum_bat.cn.bed
		"""
		# zcat {WORK}/wgs_out/{wildcards.tumor}_vs_{wildcards.normal}/battenberg/{wildcards.tumor}_battenberg_cn.vcf.gz | \
		# grep -v ^# | cut -f 1,2,8,10 | \
		# awk '{for(i=1;i<=NF;i++){if (i==3) {split($i,a,"=");printf("%s\t",a[3])} else if (i==4) {split($i,a,":");printf("%s\n", a[2])} else {printf("%s\t",$i)} }}' > \
		# {WORK}/wgs_out/{wildcards.tumor}_vs_{wildcards.normal}/tmp/norm_bat.cn.bed;
		# zcat {WORK}/wgs_out/{wildcards.tumor}_vs_{wildcards.normal}/battenberg/{wildcards.tumor}_battenberg_cn.vcf.gz | \
		# grep -v ^# | cut -f 1,2,8,11 | \
		# awk '{for(i=1;i<=NF;i++){if (i==3) {split($i,a,"=");printf("%s\t",a[3])} else if (i==4) {split($i,a,":");printf("%s\n", a[2])} else {printf("%s\t",$i)} }}' > \
		# {WORK}/wgs_out/{wildcards.tumor}_vs_{wildcards.normal}/tmp/tum_bat.cn.bed

# Block 3.5
# ASCAT_CN=/work/wgs_out/G37-B8_vs_G37-T8/WGS_G37-T8_vs_G37-B8/ascat/G37-T8.copynumber.caveman.csv
# + perl -ne '@F=(split q{,}, $_)[1,2,3,4]; $F[1]-1; print join("\t",@F)."\n";'
# + perl -ne '@F=(split q{,}, $_)[1,2,3,6]; $F[1]-1; print join("\t",@F)."\n";'



###############################################################################
#### Parallel Block 4
#### cgpPindel
#### CaVEMan
###############################################################################
rule pindel:
	input:
		normal=WORK+"/bam_out/{normal}/{normal}.bam",
		tumor=WORK+"/bam_out/{tumor}/{tumor}.bam",
		res3_1=WORK+"/wgs_out/{tumor}_vs_{normal}/timings/WGS_{tumor}_vs_{normal}.time.ascat",
		res3_2=WORK+"/wgs_out/{tumor}_vs_{normal}/timings/WGS_{tumor}_vs_{normal}.time.BRASS_cover",
		res3_3=WORK+"/wgs_out/{tumor}_vs_{normal}/timings/WGS_{tumor}_vs_{normal}.time.BRASS_input"
	output:
		timing=WORK+"/wgs_out/{tumor}_vs_{normal}/timings/WGS_{tumor}_vs_{normal}.time.cgpPindel"
	threads: 24
	shell:
		"""
		python2 ~/udocker run \
		--workdir=/  \
		--user=laran  \
		--volume={WORK}:/work  \
		cgpwgs \
		bash -c '/usr/bin/time -v pindel.pl \
		-o /work/wgs_out/{wildcards.tumor}_vs_{wildcards.normal}/WGS_{wildcards.tumor}_vs_{wildcards.normal}/pindel \
		-r /work/reference_files/genome.fa \
		-t /work/bam_out/{wildcards.tumor}/{wildcards.tumor}.bam \
		-n /work/bam_out/{wildcards.normal}/{wildcards.normal}.bam \
		-s /work/reference_files/pindel/simpleRepeats.bed.gz \
		-u /work/reference_files/pindel/pindel_np.gff3.gz \
		-f /work/reference_files/pindel/WGS_Rules.lst \
		-g /work/reference_files/vagrent/codingexon_regions.indel.bed.gz \
		-st WGS \
		-as NCBI37 \
		-sp '\''Human'\'' \
		-e NC_007605,hs37d5,GL% \
		-b /work/reference_files/pindel/HiDepth.bed.gz \
		-c {threads} \
		-sf /work/reference_files/pindel/softRules.lst \
		>& /work/wgs_out/{wildcards.tumor}_vs_{wildcards.normal}/timings/WGS_{wildcards.tumor}_vs_{wildcards.normal}.time.cgpPindel ;
		echo \''WRAPPER_EXIT: '\' $? ;
		bgzip /work/wgs_out/{wildcards.tumor}_vs_{wildcards.normal}/WGS_{wildcards.tumor}_vs_{wildcards.normal}/pindel/{wildcards.tumor}_vs_{wildcards.normal}.germline.bed;
		tabix -p bed /work/wgs_out/{wildcards.tumor}_vs_{wildcards.normal}/WGS_{wildcards.tumor}_vs_{wildcards.normal}/pindel/{wildcards.tumor}_vs_{wildcards.normal}.germline.bed.gz'
		"""

rule caveman:
	input:
		normal=WORK+"/bam_out/{normal}/{normal}.bam",
		tumor=WORK+"/bam_out/{tumor}/{tumor}.bam",
		res3_1=WORK+"/wgs_out/{tumor}_vs_{normal}/timings/WGS_{tumor}_vs_{normal}.time.ascat",
		res3_2=WORK+"/wgs_out/{tumor}_vs_{normal}/timings/WGS_{tumor}_vs_{normal}.time.BRASS_cover",
		res3_3=WORK+"/wgs_out/{tumor}_vs_{normal}/timings/WGS_{tumor}_vs_{normal}.time.BRASS_input",
		res3_4=WORK+"/wgs_out/{tumor}_vs_{normal}/tmp/tum_bat.cn.bed"
	output:
		timing=WORK+"/wgs_out/{tumor}_vs_{normal}/timings/WGS_{tumor}_vs_{normal}.time.CaVEMan"
	threads: 24
	shell:
		"""
		python2 ~/udocker run \
		--workdir=/  \
		--user=laran  \
		--volume={WORK}:/work  \
		cgpwgs \
		bash -c '/usr/bin/time -v caveman.pl \
		-r /work/reference_files/genome.fa.fai \
		-ig /work/reference_files/caveman/HiDepth.tsv \
		-b /work/reference_files/caveman/flagging \
		-ab /work/reference_files/vagrent \
		-u /work/reference_files/caveman \
		-s '\''Human'\'' \
		-sa NCBI37 \
		-t {threads} \
		-st WGS \
		-tc /work/wgs_out/{wildcards.tumor}_vs_{wildcards.normal}/tmp/tum_bat.cn.bed \
		-nc /work/wgs_out/{wildcards.tumor}_vs_{wildcards.normal}/tmp/norm_bat.cn.bed \
		-td 5 -nd 2 \
		-tb /work/bam_out/{wildcards.tumor}/{wildcards.tumor}.bam \
		-nb /work/bam_out/{wildcards.normal}/{wildcards.normal}.bam \
		-c /work/wgs_out/{wildcards.tumor}_vs_{wildcards.normal}/flag.vcf.config.WGS.ini \
		-f /work/reference_files/caveman/flagging/flag.to.vcf.convert.ini \
		-e 350000 \
		-o /work/wgs_out/{wildcards.tumor}_vs_{wildcards.normal}/WGS_{wildcards.tumor}_vs_{wildcards.normal}/caveman \
		-x NC_007605,hs37d5,GL% \
		-no-flagging \
		>& /work/wgs_out/{wildcards.tumor}_vs_{wildcards.normal}/timings/WGS_{wildcards.tumor}_vs_{wildcards.normal}.time.CaVEMan ;
		echo \''WRAPPER_EXIT: '\' $?'
		"""

###############################################################################
#### Parallel Block 5
#### BRASS
#### Pindel_annot
#### cgpFlagCaVEMan
###############################################################################
rule cgpFlagCaveman:
	input:
		normal=WORK+"/bam_out/{normal}/{normal}.bam",
		tumor=WORK+"/bam_out/{tumor}/{tumor}.bam",
		res4_1=WORK+"/wgs_out/{tumor}_vs_{normal}/timings/WGS_{tumor}_vs_{normal}.time.cgpPindel",
		res4_2=WORK+"/wgs_out/{tumor}_vs_{normal}/timings/WGS_{tumor}_vs_{normal}.time.CaVEMan",
		#res4_1_1=WORK+"/wgs_out/{tumor}_vs_{normal}/WGS_{tumor}_vs_{normal}/pindel/{tumor}_vs_{normal}.germline.bed.gz"
	output:
		timing=WORK+"/wgs_out/{tumor}_vs_{normal}/timings/WGS_{tumor}_vs_{normal}.time.cgpFlagCaVEMan"
	shell:
		"""
		python2 ~/udocker run \
		--workdir=/  \
		--user=laran  \
		--volume={WORK}:/work  \
		cgpwgs \
		bash -c '/usr/bin/time -v cgpFlagCaVEMan.pl \
		-i /work/wgs_out/{wildcards.tumor}_vs_{wildcards.normal}/WGS_{wildcards.tumor}_vs_{wildcards.normal}/caveman/{wildcards.tumor}_vs_{wildcards.normal}.muts.ids.vcf.gz \
		-o /work/wgs_out/{wildcards.tumor}_vs_{wildcards.normal}/WGS_{wildcards.tumor}_vs_{wildcards.normal}/caveman/{wildcards.tumor}_vs_{wildcards.normal}.flagged.muts.vcf \
		-s '\''Human'\'' \
		-m /work/bam_out/{wildcards.tumor}/{wildcards.tumor}.bam \
		-n /work/bam_out/{wildcards.normal}/{wildcards.normal}.bam \
		-b /work/reference_files/caveman/flagging \
		-g /work/wgs_out/{wildcards.tumor}_vs_{wildcards.normal}/WGS_{wildcards.tumor}_vs_{wildcards.normal}/pindel/{wildcards.tumor}_vs_{wildcards.normal}.germline.bed.gz \
		-umv /work/reference_files/caveman \
		-ab /work/reference_files/vagrent \
		-ref /work/reference_files/genome.fa.fai \
		-c /work/wgs_out/{wildcards.tumor}_vs_{wildcards.normal}/flag.vcf.config.WGS.ini \
		-v /work/reference_files/caveman/flagging/flag.to.vcf.convert.ini \
		>& /work/wgs_out/{wildcards.tumor}_vs_{wildcards.normal}/timings/WGS_{wildcards.tumor}_vs_{wildcards.normal}.time.cgpFlagCaVEMan ;
		echo \''WRAPPER_EXIT: '\' $?;
		bgzip /work/wgs_out/{wildcards.tumor}_vs_{wildcards.normal}/WGS_{wildcards.tumor}_vs_{wildcards.normal}/caveman/{wildcards.tumor}_vs_{wildcards.normal}.flagged.muts.vcf ;
		tabix -p vcf /work/wgs_out/{wildcards.tumor}_vs_{wildcards.normal}/WGS_{wildcards.tumor}_vs_{wildcards.normal}/caveman/{wildcards.tumor}_vs_{wildcards.normal}.flagged.muts.vcf.gz'
		"""


rule PindelAnnotate:
	input:
		res4_1=WORK+"/wgs_out/{tumor}_vs_{normal}/timings/WGS_{tumor}_vs_{normal}.time.cgpPindel",
		#res4_2=WORK+"/wgs_out/{tumor}_vs_{normal}/timings/WGS_{tumor}_vs_{normal}.time.CaVEMan"
	output:
		timing=WORK+"/wgs_out/{tumor}_vs_{normal}/timings/WGS_{tumor}_vs_{normal}.time.cgpPindel_annot"
	shell:
		"""
		python2 ~/udocker run \
		--workdir=/  \
		--user=laran  \
		--volume={WORK}:/work  \
		cgpwgs \
		bash -c '/usr/bin/time -v AnnotateVcf.pl \
		-t \
		-c /work/reference_files/vagrent/vagrent.cache.gz \
		-i /work/wgs_out/{wildcards.tumor}_vs_{wildcards.normal}/WGS_{wildcards.tumor}_vs_{wildcards.normal}/pindel/{wildcards.tumor}_vs_{wildcards.normal}.flagged.vcf.gz \
		-o /work/wgs_out/{wildcards.tumor}_vs_{wildcards.normal}/WGS_{wildcards.tumor}_vs_{wildcards.normal}/pindel/{wildcards.tumor}_vs_{wildcards.normal}.annot.vcf \
		>& /work/wgs_out/{wildcards.tumor}_vs_{wildcards.normal}/timings/WGS_{wildcards.tumor}_vs_{wildcards.normal}.time.cgpPindel_annot ;
		echo \''WRAPPER_EXIT: '\' $?'
		"""


rule Brass:
	input:
		normal=WORK+"/bam_out/{normal}/{normal}.bam",
		tumor=WORK+"/bam_out/{tumor}/{tumor}.bam",
		res3_1=WORK+"/wgs_out/{tumor}_vs_{normal}/timings/WGS_{tumor}_vs_{normal}.time.ascat",
		res3_2=WORK+"/wgs_out/{tumor}_vs_{normal}/timings/WGS_{tumor}_vs_{normal}.time.BRASS_cover",
		res3_3=WORK+"/wgs_out/{tumor}_vs_{normal}/timings/WGS_{tumor}_vs_{normal}.time.BRASS_input",
		res3_4=WORK+"/wgs_out/{tumor}_vs_{normal}/tmp/tum_bat.cn.bed",
		res4_1=WORK+"/wgs_out/{tumor}_vs_{normal}/timings/WGS_{tumor}_vs_{normal}.time.cgpPindel",
		res4_2=WORK+"/wgs_out/{tumor}_vs_{normal}/timings/WGS_{tumor}_vs_{normal}.time.CaVEMan"
	output:
		timing=WORK+"/wgs_out/{tumor}_vs_{normal}/timings/WGS_{tumor}_vs_{normal}.time.BRASS"
	shell:
		"""
		udocker run \
		--volume={WORK}:/work  \
		cgpwgs \
		bash -c '/usr/bin/time -v brass.pl \
		-j 4 -k 4 \
		-c {threads} \
		-d /work/reference_files/brass/HiDepth.bed.gz \
		-f /work/reference_files/brass/brass_np.groups.gz \
		-g /work/reference_files/genome.fa \
		-s '\''Human'\'' -as NCBI37 -pr WGS -pl ILLUMINA \
		-g_cache /work/reference_files/vagrent/vagrent.cache.gz \
		-vi /work/reference_files/brass/viral.genomic.fa.2bit \
		-mi /work/reference_files/brass/all_ncbi_bacteria \
		-b /work/reference_files/brass/500bp_windows.gc.bed.gz \
		-ct /work/reference_files/brass/CentTelo.tsv \
		-cb /work/reference_files/brass/cytoband.txt \
		-t /work/bam_out/{wildcards.tumor}/{wildcards.tumor}.bam \
		-n /work/bam_out/{wildcards.normal}/{wildcards.normal}.bam \
		-ss /work/wgs_out/{wildcards.tumor}_vs_{wildcards.normal}/battenberg/{wildcards.tumor}.samplestatistics.txt \
		-o /work/wgs_out/{wildcards.tumor}_vs_{wildcards.normal}/WGS_{wildcards.tumor}_vs_{wildcards.normal}/brass \
		>& /work/wgs_out/{wildcards.tumor}_vs_{wildcards.normal}/timings/WGS_{wildcards.tumor}_vs_{wildcards.normal}.time.BRASS ;
		echo \''WRAPPER_EXIT: '\' $?'
		"""
		# -ss /work/wgs_out/{wildcards.tumor}_vs_{wildcards.normal}/WGS_{wildcards.tumor}_vs_{wildcards.normal}/ascat/{wildcards.tumor}.samplestatistics.txt \

###############################################################################
#### Parallel Block 6
#### CaVEMan_annot
#### VerifyBam tumor
###############################################################################
rule VerifyBam_tumor:
	input:
		tumor=WORK+"/bam_out/{tumor}/{tumor}.bam",
		res3_1=WORK+"/wgs_out/{tumor}_vs_{normal}/timings/WGS_{tumor}_vs_{normal}.time.ascat",
		res3_4=WORK+"/wgs_out/{tumor}_vs_{normal}/tmp/tum_bat.cn.bed",
		res1_4=WORK+"/wgs_out/{tumor}_vs_{normal}/timings/WGS_{tumor}_vs_{normal}.time.verify_WT"
	output:
		timing=WORK+"/wgs_out/{tumor}_vs_{normal}/timings/WGS_{tumor}_vs_{normal}.time.verify_MT"
	threads: 24
	shell:
		"""
		python2 ~/udocker run \
		--workdir=/  \
		--user=laran  \
		--volume={WORK}:/work  \
		cgpwgs \
		bash -c '/usr/bin/time -v verifyBamHomChk.pl \
		-d 25 \
		-o /work/wgs_out/{wildcards.tumor}_vs_{wildcards.normal}/WGS_{wildcards.tumor}/contamination \
		-b /work/bam_out/{wildcards.tumor}/{wildcards.tumor}.bam \
		-t {threads} \
		-a /work/wgs_out/{wildcards.tumor}_vs_{wildcards.normal}/battenberg/{wildcards.tumor}_battenberg_copynumber.csv \
		-j /work/wgs_out/{wildcards.tumor}_vs_{wildcards.normal}/WGS_{wildcards.tumor}/contamination/result.json \
		-s /work/reference_files/verifyBamID_snps.vcf.gz \
		>& /work/wgs_out/{wildcards.tumor}_vs_{wildcards.normal}/timings/WGS_{wildcards.tumor}_vs_{wildcards.normal}.time.verify_MT ;
		echo \''WRAPPER_EXIT: '\' $?'
		"""
		#-a /work/wgs_out/{wildcards.tumor}_vs_{wildcards.normal}/WGS_{wildcards.tumor}_vs_{wildcards.normal}/ascat/{wildcards.tumor}.copynumber.caveman.csv \

rule caveman_annot:
	input:
		res5_1=WORK+"/wgs_out/{tumor}_vs_{normal}/timings/WGS_{tumor}_vs_{normal}.time.cgpFlagCaVEMan",
		res5_2=WORK+"/wgs_out/{tumor}_vs_{normal}/timings/WGS_{tumor}_vs_{normal}.time.cgpPindel_annot",
		#res5_3=WORK+"/wgs_out/{tumor}_vs_{normal}/timings/WGS_{tumor}_vs_{normal}.time.BRASS",
		#res5_1_1=WORK+"/wgs_out/{tumor}_vs_{normal}/WGS_{tumor}_vs_{normal}/caveman/{tumor}_vs_{normal}.flagged.muts.vcf.gz"
	output:
		timing=WORK+"/wgs_out/{tumor}_vs_{normal}/timings/WGS_{tumor}_vs_{normal}.time.CaVEMan_annot"
	shell:
		"""
		python2 ~/udocker run \
		--workdir=/  \
		--user=laran  \
		--volume={WORK}:/work  \
		cgpwgs \
		bash -c '/usr/bin/time -v AnnotateVcf.pl \
		-t \
		-c /work/reference_files/vagrent/vagrent.cache.gz \
		-i /work/wgs_out/{wildcards.tumor}_vs_{wildcards.normal}/WGS_{wildcards.tumor}_vs_{wildcards.normal}/caveman/{wildcards.tumor}_vs_{wildcards.normal}.flagged.muts.vcf.gz \
		-o /work/wgs_out/{wildcards.tumor}_vs_{wildcards.normal}/WGS_{wildcards.tumor}_vs_{wildcards.normal}/caveman/{wildcards.tumor}_vs_{wildcards.normal}.annot.muts.vcf \
		>& /work/wgs_out/{wildcards.tumor}_vs_{wildcards.normal}/timings/WGS_{wildcards.tumor}_vs_{wildcards.normal}.time.CaVEMan_annot ;
		echo \''WRAPPER_EXIT: '\' $?'
		"""

rule caveman_annot_germline:
	input:
		# res5_1=WORK+"/wgs_out/{tumor}_vs_{normal}/timings/WGS_{tumor}_vs_{normal}.time.cgpFlagCaVEMan",
		# res5_2=WORK+"/wgs_out/{tumor}_vs_{normal}/timings/WGS_{tumor}_vs_{normal}.time.cgpPindel_annot",
		# res5_3=WORK+"/wgs_out/{tumor}_vs_{normal}/timings/WGS_{tumor}_vs_{normal}.time.BRASS",
		res4_2=WORK+"/wgs_out/{tumor}_vs_{normal}/timings/WGS_{tumor}_vs_{normal}.time.CaVEMan",
	output:
		timing=WORK+"/wgs_out/{tumor}_vs_{normal}/timings/WGS_{tumor}_vs_{normal}.time.CaVEMan_annot_germline"
	shell:
		"""
		python2 ~/udocker run \
		--workdir=/  \
		--user=laran  \
		--volume={WORK}:/work  \
		cgpwgs \
		bash -c '/usr/bin/time -v AnnotateVcf.pl \
		-t \
		-c /work/reference_files/vagrent/vagrent.cache.gz \
		-i /work/wgs_out/{wildcards.tumor}_vs_{wildcards.normal}/WGS_{wildcards.tumor}_vs_{wildcards.normal}/caveman/{wildcards.tumor}_vs_{wildcards.normal}.snps.ids.vcf.gz \
		-o /work/wgs_out/{wildcards.tumor}_vs_{wildcards.normal}/WGS_{wildcards.tumor}_vs_{wildcards.normal}/caveman/{wildcards.tumor}_vs_{wildcards.normal}.snps.ids.annot.vcf \
		>& /work/wgs_out/{wildcards.tumor}_vs_{wildcards.normal}/timings/WGS_{wildcards.tumor}_vs_{wildcards.normal}.time.CaVEMan_annot_germline ;
		echo \''WRAPPER_EXIT: '\' $?'
		"""



###############################################################################
#### Execution of battenberg copynumber (Manual)
###############################################################################
rule battenberg:
	input:
		normal=WORK+"/bam_out/{normal}/{normal}.bam",
		tumor=WORK+"/bam_out/{tumor}/{tumor}.bam",
		res2_1=WORK+"/wgs_out/{tumor}_vs_{normal}/timings/WGS_{tumor}_vs_{normal}.time.CaVEMan_split",
		res2_2=WORK+"/wgs_out/{tumor}_vs_{normal}/timings/WGS_{tumor}_vs_{normal}.time.alleleCount"
	output:
		timing=WORK+"/wgs_out/{tumor}_vs_{normal}/timings/WGS_{tumor}_vs_{normal}.time.battenberg"
	threads: 24
	shell:
		"""
		python2 ~/udocker run \
		--workdir=/  \
		--user=laran  \
		--volume={WORK}:/work  \
		cgpwgs201 \
		bash -c '/usr/bin/time -v battenberg.pl \
		-o /work/wgs_out/{wildcards.tumor}_vs_{wildcards.normal}/battenberg \
		-u /work/reference_files/battenberg/1000genomesloci \
		-e /work/reference_files/battenberg/impute/impute_info.txt \
		-c /work/reference_files/battenberg/probloci.txt \
		-r /work/reference_files/genome.fa.fai \
		-ig /work/reference_files/battenberg/ignore_contigs.txt \
		-gc /work/reference_files/battenberg/battenberg_wgs_gc_correction_1000g \
		-ge XX \
		-tb /work/bam_out/{wildcards.tumor}/{wildcards.tumor}.bam \
		-nb /work/bam_out/{wildcards.normal}/{wildcards.normal}.bam \
		-rs 'Human' \
		-ra NCBI37 \
		-nl 50 \
		-t {threads} \
		>& /work/wgs_out/{wildcards.tumor}_vs_{wildcards.normal}/timings/WGS_{wildcards.tumor}_vs_{wildcards.normal}.time.battenberg ;
		echo '\''WRAPPER_EXIT: '\''$?';
		"""
		# cp -r {WORK}/wgs_out/{wildcards.tumor}_vs_{wildcards.normal}/WGS_{wildcards.tumor}_vs_{wildcards.normal}/battenberg {WORK}/wgs_out/{wildcards.tumor}_vs_{wildcards.normal}/battenberg; \
		#GRCh37d5
		#cp -r /work/wgs_out/{wildcards.tumor}_vs_{wildcards.normal}/WGS_{wildcards.tumor}_vs_{wildcards.normal}/battenberg ../battenberg_bak; \
		# zcat /work/wgs_out/{{wildcards.tumor}}_vs_{{wildcards.normal}}/battenberg/{{wildcards.tumor}}_battenberg_cn.vcf.gz | \
		# grep -v ^# | cut -f 1,2,8,10 | \
		# awk '{for(i=1;i<=NF;i++){if (i==3) {split($i,a,"=");printf("%s\t",a[3])} else if (i==4) {split($i,a,":");printf("%s\n", a[2])} else {printf("%s\t",$i)} }}' > \
		# /work/wgs_out/{wildcards.tumor}_vs_{wildcards.normal}/tmp/norm_bat.cn.bed;
		# zcat /work/wgs_out/{wildcards.tumor}_vs_{wildcards.normal}/battenberg/{wildcards.tumor}_battenberg_cn.vcf.gz | \
		# grep -v ^# | cut -f 1,2,8,11 | \
		# awk '{for(i=1;i<=NF;i++){if (i==3) {split($i,a,"=");printf("%s\t",a[3])} else if (i==4) {split($i,a,":");printf("%s\n", a[2])} else {printf("%s\t",$i)} }}' > \
		# /work/wgs_out/{wildcards.tumor}_vs_{wildcards.normal}/tmp/tum_bat.cn.bed
