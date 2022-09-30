__title__ = "Pipeline for Somatic Variant Calling with Mutect2 - Kristina's project"
__author__ = "Lars Andersen <larsmew@gmail.com>"
__date__ = "10/09/2019"
__version__ = "1.0"

import time, os, sys, glob

#########################################################
####                       Input                     ####
#########################################################
# Matched tumour-normal samples information
configfile: "somatic_matched_samples.yaml"
# configfile: "somatic_matched_samples_oneSample.yaml"

# Explicit paths for external input files
ref = "/work/sdukoldby/resources/hg38/Homo_sapiens_assembly38.fasta"
# gnomead = "/work/sduvarcall/knownSNPs/gnomead/af-only-gnomad.raw.sites.b37.vcf.gz"
gnomad = "/work/sdukoldby/resources/hg38/af-only-gnomad.hg38.vcf.gz"
interval_list = "/work/sdukoldby/resources/hg38/MedExome_hg38_capture_targets.interval_list"
target_regions = "/work/sdukoldby/resources/hg38/target_regions/"
common_variants = "/work/sdukoldby/resources/hg38/small_exac_common_3.hg38.vcf.gz"
# Panel of normals location
pon_location = "somatic_panel_of_normals/"

# Sample information
PAIR, = glob_wildcards("{sample}_tumor_tagseq-medexome.connor.recalibrated.bam")
# PAIR = ["ECV2-4-plasma180205"]
PAIR = [pair+"_tumor" for pair in PAIR]
print(PAIR)

NORMALS, = glob_wildcards(pon_location+"{normal}.bam")
#NORMALS = ["ECV2-35-blod_normal_tagseq-medexome.recalibrated"]
print(NORMALS)
# sys.exit()

INTERVALS, = glob_wildcards(target_regions+"{interval}.interval_list")
INTERVALS = sorted(INTERVALS)
print(INTERVALS)
# sys.exit()

# SAMPLES = "ECV2-29-blod_normal_tagseq-medexome.connor.recalibrated_for_pon.vcf"

#########################################################
####                      Output                     ####
#########################################################
log_file = "log_file_somatic.txt"
output_somatic = "somatic_variants/"


#########################################################
####                       Setup                     ####
#########################################################
# Timing
totim = time.time()
timeFormat = "%Y_%m_%d:%X" # year, month, day, time H:M:S

# Memory
# mem = "-Xmx12g" # login nodes - be careful not running too many jobs at once!
# mem = "-Xmx24g" # slim nodes
mem = "-Xmx32g" # Fat nodes
# mem = "-Xmx64g"

#########################################################
####  Define workflow start, stop and error actions  ####
#########################################################
onstart:
	# shell("echo $(head /work/sduvarcall/exome_analysis/pipeline_version.txt -n 1) >> {log_file}")
	# shell("echo 'Started execution of pipeline:' $(date +'%Y-%m-%d %H:%M:%S') >> {log_file}")
	shell("mkdir -p "+output_somatic)
	shell("mkdir -p somatic_panel_of_normals/split")
	shell("mkdir -p "+output_somatic+"split/")
	shell("mkdir -p split_f1r2/")


#########################################################
####                  Run All Rules                  ####
#########################################################
'''
Rule all
'''
print([expand(output_somatic+"{tumor}_vs_{normal}_tagseq-medexome_somatic_mutect2_filtered.vcf",
			normal=config[pair]["normal"],
			tumor=config[pair]["tumor"]) for pair in PAIR])
# sys.exit()

rule all_pairs:
	input:
		[expand(output_somatic+"{tumor}_vs_{normal}_tagseq-medexome_somatic_mutect2_filtered.vcf",
			normal=config[fam]["normal"],
			tumor=config[fam]["tumor"]) for fam in PAIR]
		# pon_location+"pon.vcf.gz"
		# pon_location+"ECV2-29-blod_normal_tagseq-medexome.connor.recalibrated_for_pon.vcf"
		# pon_location+"ECV2-31-blod_normal_tagseq-medexome.connor.recalibrated_for_pon.vcf"
		# pon_location+"ECV2-35-blod_normal_tagseq-medexome.connor.recalibrated_for_pon.vcf"
		# pon_location+"ECV2-4-blod_normal_tagseq-medexome.connor.recalibrated_for_pon.vcf"
		# pon_location+"ECV2-8-blod_normal_tagseq-medexome.connor.recalibrated_for_pon.vcf"
		# pon_location+"PC1-10-blod_normal_tagseq-medexome.connor.recalibrated_for_pon.vcf"
		# pon_location+"PC1-14-blod_normal_tagseq-medexome.connor.recalibrated_for_pon.vcf"
		# pon_location+"PC1-18-blod_normal_tagseq-medexome.connor.recalibrated_for_pon.vcf"
		# expand(pon_location+"{sample}_for_pon.vcf", sample=NORMALS)
		# [expand("{tumor}_vs_{normal}_tagseq-medexome_somatic_mutect2.vcf",
		# 			normal=config[fam]["normal"],
		# 			tumor=config[fam]["tumor"]) for fam in PAIR]


#########################################################
####       Create Somatic Panel of Normals           ####
#########################################################
'''
Run Mutect2 in tumor-only mode for each normal sample
'''
rule Mutect2_tumor_only_pon:
	input:
		bam=pon_location+"{sample}.bam",
		intervals=target_regions+"{interval}.interval_list"
	output:
		vcf=pon_location+"split/{sample}_for_pon__{interval}__split.vcf",
		idx=pon_location+"split/{sample}_for_pon__{interval}__split.vcf.idx"
	threads: 24
	shell:
		"""
		gatk --java-options {mem} Mutect2 \
		-R {ref} \
		-I {input.bam} \
		-max-mnp-distance 0 \
		--native-pair-hmm-threads {threads} \
		-L {input.intervals} \
		-O {output.vcf}
		"""

rule merge_normal_vcf:
	input:
		vcf_subfile=expand(pon_location+"split/{{sample}}_for_pon__{interval}__split.vcf", interval=INTERVALS)
	output:
		vcf=pon_location+"{sample}_for_pon.vcf",
		idx=pon_location+"{sample}_for_pon.vcf.idx"
	params:
		vcf_subfile=expand("-I "+pon_location+"split/{{sample}}_for_pon__{interval}__split.vcf", interval=INTERVALS)
	shell:
		"""
		gatk --java-options {mem} GatherVcfs \
		{params.vcf_subfile} \
		-O {output.vcf}
		"""

rule GenomicsDB:
	input:
		vcf=expand(pon_location+"{sample}_for_pon.vcf.gz", sample=NORMALS)
	output:
		db=directory(pon_location+"pon_db")
	params:
		vcf=expand("-V "+pon_location+"{sample}_for_pon.vcf.gz", sample=NORMALS)
	shell:
		"""
		gatk --java-options {mem} GenomicsDBImport \
			-R {ref} \
			-L {interval_list} \
			--genomicsdb-workspace-path {output.db} \
			--max-num-intervals-to-import-in-parallel {threads} \
			--merge-input-intervals true \
			{params.vcf}
		"""

'''
Combine the normal calls using CreateSomaticPanelOfNormals.
'''
rule CreateSomaticPanelOfNormals:
	input:
		db=pon_location+"pon_db"
		# vcf=expand(pon_location+"{sample}_for_pon.vcf", sample=NORMALS)
	output:
		pon=pon_location+"pon.vcf.gz"
	# params:
	# 	vcf=expand("-V "+pon_location+"{sample}_for_pon.vcf", sample=NORMALS)
	shell:
		"""
		gatk --java-options {mem} CreateSomaticPanelOfNormals \
		-R {ref} \
		-V gendb://{input.db} \
		-O {output}
		"""
		# {params.vcf} \


##########################################################
####  Call Somatic Variants using Mutect2 on matched  ####
####   Tumor-Normal samples on per chromesome basis   ####
##########################################################
'''
Mutect2 on matched Tumor-Normal samples
'''
rule Mutect2_matched:
	input:
		normal="{normal}_tagseq-medexome.connor.recalibrated.bam",
		tumor="{tumor}_tagseq-medexome.connor.recalibrated.bam",
		pon=pon_location+"pon.vcf.gz",
		intervals=target_regions+"{interval}.interval_list"
	output:
		vcf=output_somatic+"split/{tumor}_vs_{normal}_tagseq-medexome_somatic_mutect2__{interval}__split.vcf",
		idx=output_somatic+"split/{tumor}_vs_{normal}_tagseq-medexome_somatic_mutect2__{interval}__split.vcf.idx",
		vcf_stats=output_somatic+"split/{tumor}_vs_{normal}_tagseq-medexome_somatic_mutect2__{interval}__split.vcf.stats",
		f1r2="split_f1r2/{tumor}_vs_{normal}_f1r2__{interval}__split.tar.gz"
	threads: 24
	shell:
		"""
		gatk --java-options {mem} Mutect2 \
		-R {ref} \
		-I {input.tumor} \
		-I {input.normal} \
        -normal {wildcards.normal}_tagseq-medexome \
		-pon {input.pon} \
		--germline-resource {gnomad} \
		--af-of-alleles-not-in-resource 0.0000025 \
		--disable-read-filter MateOnSameContigOrNoMappedMateReadFilter \
		--native-pair-hmm-threads {threads} \
		--f1r2-tar-gz {output.f1r2} \
        -bamout bam_split/{wildcards.tumor}__{wildcards.interval}_mutect2.bam \
		-L {input.intervals} \
		-O {output.vcf}
		"""
		# -bamout bam_split/{wildcards.tumor}__{wildcards.interval}_mutect2.bam
		# -tumor {wildcards.tumor} \
		# -normal {wildcards.normal} \

rule merge_somatic_vcf:
	input:
		vcf_subfile=expand(output_somatic+"split/{{tumor}}_vs_{{normal}}_tagseq-medexome_somatic_mutect2__{interval}__split.vcf", interval=INTERVALS)
	output:
		vcf=output_somatic+"{tumor}_vs_{normal}_tagseq-medexome_somatic_mutect2.vcf",
		idx=output_somatic+"{tumor}_vs_{normal}_tagseq-medexome_somatic_mutect2.vcf.idx",
	params:
		vcf_subfile=expand("-I "+output_somatic+"split/{{tumor}}_vs_{{normal}}_tagseq-medexome_somatic_mutect2__{interval}__split.vcf", interval=INTERVALS)
	shell:
		"""
		gatk --java-options {mem} GatherVcfs \
		{params.vcf_subfile} \
		-O {output.vcf}
		"""
	
rule merge_somatic_vcf_stats:
	input:
		vcf_subfile=expand(output_somatic+"split/{{tumor}}_vs_{{normal}}_tagseq-medexome_somatic_mutect2__{interval}__split.vcf.stats", interval=INTERVALS)
	output:
		vcf_stats=output_somatic+"{tumor}_vs_{normal}_tagseq-medexome_somatic_mutect2.vcf.stats",
	params:
		vcf_subfile=expand("-stats "+output_somatic+"split/{{tumor}}_vs_{{normal}}_tagseq-medexome_somatic_mutect2__{interval}__split.vcf.stats", interval=INTERVALS)
	shell:
		"""
		gatk --java-options {mem} MergeMutectStats \
		{params.vcf_subfile} \
		-O {output.vcf_stats}
		"""

#########################################################
####          Learn Read Orientation Bias            ####
#########################################################
rule learnReadOrientationModel:
	input:
		f1r2=expand("split_f1r2/{{tumor}}_vs_{{normal}}_f1r2__{interval}__split.tar.gz", interval=INTERVALS)
	output:
		f1r2_model=output_somatic+"{tumor}_vs_{normal}_read-orientation-model.tar.gz"
	params:
		f1r2=expand("-I split_f1r2/{{tumor}}_vs_{{normal}}_f1r2__{interval}__split.tar.gz", interval=INTERVALS)
	shell:
		"""
		gatk LearnReadOrientationModel \
		{params.f1r2} \
		-O {output}
		"""
		# -I {input} \

#########################################################
####           Create Contamination table            ####
#########################################################
### Sæt sammen til en metode, som tager højde for med og uden normal... (Se gatk4 scripts)
rule GetPileupSummaries_normal:
	input:
		bam="{normal}_tagseq-medexome.connor.recalibrated.bam"
	output:
		pileup=output_somatic+"{normal}_normal_pileup.table"
	shell:
		"""
		gatk --java-options {mem} GetPileupSummaries \
		-I {input.bam} \
		-V {common_variants} \
		-L {common_variants} \
		-O {output}
		"""

rule GetPileupSummaries_tumor:
	input:
		bam="{tumor}_tagseq-medexome.connor.recalibrated.bam"
	output:
		pileup=output_somatic+"{tumor}_tumor_pileup.table"
	shell:
		"""
		gatk --java-options {mem} GetPileupSummaries \
		-I {input.bam} \
		-V {common_variants} \
		-L {common_variants} \
		-O {output}
		"""


rule CalculateContamination:
	input:
		normal=output_somatic+"{normal}_normal_pileup.table",
		tumor=output_somatic+"{tumor}_tumor_pileup.table"
	output:
		contamination=output_somatic+"{tumor}_vs_{normal}_contamination.table"
	shell:
		"""
		gatk --java-options {mem} CalculateContamination \
		-I {input.tumor} \
		-matched {input.normal} \
		-O {output}
		"""


#########################################################
####              Filter Mutect2 Calls               ####
#########################################################
rule FilterMutectCalls:
	input:
		vcf=output_somatic+"{tumor}_vs_{normal}_tagseq-medexome_somatic_mutect2.vcf",
		idx=output_somatic+"{tumor}_vs_{normal}_tagseq-medexome_somatic_mutect2.vcf.idx",
		vcf_stats=output_somatic+"{tumor}_vs_{normal}_tagseq-medexome_somatic_mutect2.vcf.stats",
		contamination=output_somatic+"{tumor}_vs_{normal}_contamination.table",
		read_orientation=output_somatic+"{tumor}_vs_{normal}_read-orientation-model.tar.gz"
	output:
		vcf=output_somatic+"{tumor}_vs_{normal}_tagseq-medexome_somatic_mutect2_filtered.vcf",
		idx=output_somatic+"{tumor}_vs_{normal}_tagseq-medexome_somatic_mutect2_filtered.vcf.idx"
	shell:
		"""
		gatk --java-options {mem} FilterMutectCalls \
		-R {ref} \
		-V {input.vcf} \
		--contamination-table {input.contamination} \
		--orientation-bias-artifact-priors {input.read_orientation} \
		-L {interval_list} \
		-O {output.vcf}
		"""
