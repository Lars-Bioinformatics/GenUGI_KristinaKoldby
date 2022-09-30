__title__ = "Pipeline for Somatic JOINT Variant Calling with Mutect2 - Kristina's project"
__author__ = "Lars Andersen <larsmew@gmail.com> and Kristina Koldby"
__date__ = "16/12/2019"
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
# gnomad = "/work/sduvarcall/knownSNPs/gnomead/af-only-gnomad.raw.sites.b37.vcf.gz"
gnomad = "/work/sdukoldby/resources/hg38/af-only-gnomad.hg38.vcf.gz"
interval_list = "/work/sdukoldby/resources/hg38/MedExome_hg38_capture_targets.interval_list"
target_regions = "/work/sdukoldby/resources/hg38/target_regions/"
common_variants = "/work/sdukoldby/resources/hg38/small_exac_common_3.hg38.vcf.gz"
# Panel of normals location
# pon_location = "somatic_panel_of_normals/"
pon_location = "/work/sdukoldby/data/G45-2016_genugi/exome_fastq/exome_fastq_merged/hg38/Connor/somatic_panel_of_normals/"

# SAMPLES = ["ECV2-4"]
# TUMORS = ["ECV2-4-biopsi-merged_tumor_tagseq-medexome-deep-seq","ECV2-4-op-merged_tumor_tagseq-medexome-deep-seq","ECV2-4-op-1401_tumor_tagseq-medexome-deep-seq","ECV2-4-op18-merged_tumor_tagseq-medexome-deep-seq","ECV2-4-plasma170503_tumor_tagseq-medexome","ECV2-4-plasma170712_tumor_tagseq-medexome","ECV2-4-plasma180205_tumor_tagseq-medexome"]
# NORMALS = ["ECV2-4-blod_normal_tagseq-medexome"]
# SAMPLES = ["ECV2-8"]
# TUMORS = ["ECV2-8-biopsi-merged_tumor_tagseq-medexome-deep-seq","ECV2-8-op-merged_tumor_tagseq-medexome-deep-seq","ECV2-8-plasma170522_tumor_tagseq-medexome","ECV2-8-plasma170728_tumor_tagseq-medexome","ECV2-8-plasma180806_tumor_tagseq-medexome"]
# NORMALS = ["ECV2-8-blod_normal_tagseq-medexome"]
SAMPLES = ["ECV2-29"]
TUMORS = ["ECV2-29-biopsi-merged_tumor_tagseq-medexome-deep-seq","ECV2-29-opA-merged_tumor_tagseq-medexome-deep-seq","ECV2-29-opB-merged_tumor_tagseq-medexome-deep-seq","ECV2-29-plasma171124_tumor_tagseq-medexome","ECV2-29-plasma180119_tumor_tagseq-medexome","ECV2-29-plasma180619_tumor_tagseq-medexome"]
NORMALS = ["ECV2-29-blod_normal_tagseq-medexome"]
# SAMPLES = ["ECV2-31"]
# TUMORS = ["ECV2-31-biopsi-merged_tumor_tagseq-medexome-deep-seq","ECV2-31-opD-merged_tumor_tagseq-medexome-deep-seq","ECV2-31-opE-merged_tumor_tagseq-medexome-deep-seq","ECV2-31-plasma171215_tumor_tagseq-medexome","ECV2-31-plasma180202_tumor_tagseq-medexome"]
# NORMALS = ["ECV2-31-blod_normal_tagseq-medexome"]
# SAMPLES = ["ECV2-35"]
# TUMORS = ["ECV2-35-biopsi-merged_tumor_tagseq-medexome-deep-seq","ECV2-35-op1to4_tumor_tagseq-medexome-deep-seq","ECV2-35-op6to7_tumor_tagseq-medexome-deep-seq","ECV2-35-plasma180102_tumor_tagseq-medexome","ECV2-35-plasma180316_tumor_tagseq-medexome","ECV2-35-plasma180601_tumor_tagseq-medexome"]
# NORMALS = ["ECV2-35-blod_normal_tagseq-medexome"]
# SAMPLES = ["PC1-10"]
# TUMORS = ["PC1-10-op-merged_tumor_tagseq-medexome-deep-seq","PC1-10-plasma170928_tumor_tagseq-medexome"]
# NORMALS = ["PC1-10-blod_normal_tagseq-medexome"]
# SAMPLES = ["PC1-14"]
# TUMORS = ["PC1-14-op-merged_tumor_tagseq-medexome-deep-seq","PC1-14-plasma180111_tumor_tagseq-medexome","PC1-14-plasma180907_tumor_tagseq-medexome"]
# NORMALS = ["PC1-14-blod_normal_tagseq-medexome"]
# SAMPLES = ["PC1-18"]
# TUMORS = ["PC1-18-op-merged_tumor_tagseq-medexome-deep-seq","PC1-18-plasma180219_tumor_tagseq-medexome"]
# NORMALS = ["PC1-18-blod_normal_tagseq-medexome"]


print(SAMPLES)
print(TUMORS)
print(NORMALS)

# Sample information
# PAIR, = glob_wildcards("{sample}_tumor.connor.recalibrated.bam")
# # PAIR = ["ECV2-4-plasma180205"]
# PAIR = [pair+"_tumor" for pair in PAIR]
# print(PAIR)

INTERVALS, = glob_wildcards(target_regions+"{interval}.interval_list")
INTERVALS = sorted(INTERVALS)
print(INTERVALS)
# sys.exit()

# SAMPLES = "ECV2-29-blod_normal.connor.recalibrated_for_pon.vcf"

#########################################################
####                      Output                     ####
#########################################################
log_file = "log_file_somatic.txt"
output_somatic = "mutect2_somatic_joint_calling/"


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
    shell("mkdir -p "+output_somatic)
    shell("mkdir -p "+output_somatic+"split/")
    shell("mkdir -p "+output_somatic+"split_f1r2/")


#########################################################
####                  Run All Rules                  ####
#########################################################
'''
Rule all
'''
print(expand(output_somatic+"{sample}_somatic_mutect2_filtered.vcf", sample=SAMPLES))
# sys.exit()

rule all_pairs:
    input:
        expand(output_somatic+"{sample}_somatic_mutect2_filtered.vcf", sample=SAMPLES),
        # expand(output_somatic+"{sample}_somatic_mutect2.bam", sample=SAMPLES)
        # expand(output_somatic+"{sample}_somatic_mutect2.vcf", sample=SAMPLES),
        # expand(output_somatic+"{sample}_merged_contamination.table", sample=SAMPLES)


##########################################################
####  Call Somatic Variants using Mutect2 on matched  ####
####   Tumor-Normal samples on per chromesome basis   ####
##########################################################
'''
Mutect2 on matched Tumor-Normal samples
'''
rule Mutect2_matched:
    input:
        # normal=config["{sample}"]["normal"]+".connor.recalibrated.bam",
        # tumor="{tumor}.connor.recalibrated.bam",
        normal=expand("{normal}.connor.recalibrated.bam", normal=NORMALS),
        tumors=expand("{tumor}.connor.recalibrated.bam", tumor=TUMORS),
        pon=pon_location+"pon.vcf.gz",
        intervals=target_regions+"{interval}.interval_list"
    output:
        vcf=output_somatic+"split/{sample}_somatic_mutect2__{interval}__split.vcf",
        idx=output_somatic+"split/{sample}_somatic_mutect2__{interval}__split.vcf.idx",
        vcf_stats=output_somatic+"split/{sample}_somatic_mutect2__{interval}__split.vcf.stats",
        # bam_subfile="bam_split/{sample}__{interval}_mutect2.bam",
        f1r2=output_somatic+"split_f1r2/{sample}_f1r2__{interval}__split.tar.gz"
    params:
        tumors=expand("-I {tumor}.connor.recalibrated.bam", tumor=TUMORS),
        normal_name=expand("-normal {normal}", normal=NORMALS)
    threads: 24
    shell:
        """
        gatk --java-options {mem} Mutect2 \
        -R {ref} \
        {params.tumors} \
        -I {input.normal} \
        {params.normal_name} \
        -pon {input.pon} \
        --germline-resource {gnomad} \
        --af-of-alleles-not-in-resource 0.0000025 \
        --disable-read-filter MateOnSameContigOrNoMappedMateReadFilter \
        --native-pair-hmm-threads {threads} \
        --f1r2-tar-gz {output.f1r2} \
        -L {input.intervals} \
        -O {output.vcf}
        """
        # -bamout bam_split/{wildcards.tumor}__{wildcards.interval}_mutect2.bam
        # -tumor {wildcards.tumor} \
        # -normal {wildcards.normal} \
        # -bamout {output.bam_subfile} \
        # -I {input.tumor} \
        # -normal {wildcards.normal} \
        # -bamout {output.bam_subfile} \ ## Debugging

rule merge_somatic_vcf:
    input:
        vcf_subfile=expand(output_somatic+"split/{{sample}}_somatic_mutect2__{interval}__split.vcf", interval=INTERVALS)
    output:
        vcf=output_somatic+"{sample}_somatic_mutect2.vcf",
        idx=output_somatic+"{sample}_somatic_mutect2.vcf.idx",
    params:
        vcf_subfile=expand("-I "+output_somatic+"split/{{sample}}_somatic_mutect2__{interval}__split.vcf", interval=INTERVALS)
    shell:
        """
        gatk --java-options {mem} GatherVcfs \
        {params.vcf_subfile} \
        -O {output.vcf}
        """

## FOR DEBUGGING
rule merge_mutect2_bam:
  input:
        bam_subfile=expand("bam_split/{{sample}}__{interval}_mutect2.bam", interval=INTERVALS)
  output:
        bam=output_somatic+"{sample}_somatic_mutect2.bam"
  params:
        bam_subfile=expand("-I bam_split/{{sample}}__{interval}_mutect2.bam", interval=INTERVALS)
  shell:
        """
        gatk --java-options {mem} MergeSamFiles \
        {params.bam_subfile} \
        --VALIDATION_STRINGENCY LENIENT \
        --USE_THREADING true \
        -O {output.bam} \
        --CREATE_INDEX true
        """

rule merge_somatic_vcf_stats:
    input:
        vcf_subfile=expand(output_somatic+"split/{{sample}}_somatic_mutect2__{interval}__split.vcf.stats", interval=INTERVALS)
    output:
        vcf_stats=output_somatic+"{sample}_somatic_mutect2.vcf.stats",
    params:
        vcf_subfile=expand("-stats "+output_somatic+"split/{{sample}}_somatic_mutect2__{interval}__split.vcf.stats", interval=INTERVALS)
    shell:
        """
        gatk --java-options {mem} MergeMutectStats \
        {params.vcf_subfile} \
        -O {output.vcf_stats}
        """

#########################################################
####          Learn Read Orientation Bias            ####
####                                                 ####
#### Note: Used to fix orientation bias artifacts    ####
####       from formalin Formalin-Fixed Paraffin-    ####
####       Embedded (FFPE) samples - i.e. not needed ####
####       for frozen tissue                         ####
#########################################################
rule learnReadOrientationModel:
    input:
        f1r2=expand(output_somatic+"split_f1r2/{{sample}}_f1r2__{interval}__split.tar.gz", interval=INTERVALS)
    output:
        f1r2_model=output_somatic+"{sample}_read-orientation-model.tar.gz"
    params:
        f1r2=expand("-I "+output_somatic+"split_f1r2/{{sample}}_f1r2__{interval}__split.tar.gz", interval=INTERVALS)
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
rule GetPileupSummaries:
    input:
        bam="{sample}.connor.recalibrated.bam"
    output:
        pileup=output_somatic+"{sample}_pileup.table"
    shell:
        """
        gatk --java-options {mem} GetPileupSummaries \
        -I {input.bam} \
        -V {common_variants} \
        -L {common_variants} \
        -O {output}
        """

# rule GetPileupSummaries_tumor:
#     input:
#         bam="{tumor}.connor.recalibrated.bam"
#     output:
#         pileup=output_somatic+"{tumor}_pileup.table"
#     shell:
#         ""
#         gatk --java-options {mem} GetPileupSummaries \
#         -I {input.bam} \
#         -V {common_variants} \
#         -L {common_variants} \
#         -O {output}
#         ""


rule CalculateContamination:
    input:
        normal=output_somatic+"{normal}_pileup.table",
        tumor=output_somatic+"{tumor}_pileup.table"
    output:
        contamination=output_somatic+"{tumor}_vs_{normal}_contamination.table"
    shell:
        """
        gatk --java-options {mem} CalculateContamination \
        -I {input.tumor} \
        -matched {input.normal} \
        -O {output}
        """

rule JoinContaminationTables:
    input:
        cont_tables=expand(output_somatic+"{tumor}_vs_{normal}_contamination.table", normal=NORMALS, tumor=TUMORS)
    output:
        joined_table=output_somatic+"{sample}_merged_contamination.table"
    params:
        sample_name=expand('{sample}', sample=SAMPLES)
    shell:
        """
        echo -e 'sample\tcontamination\terror' > {output};
        awk -F'\t' 'FNR == 2' {input} >> {output}
        """
        # awk -F'\t' 'FNR == 2' {params.sample_name}*contamination.table >> {output}


#########################################################
####              Filter Mutect2 Calls               ####
#########################################################
# rule FilterMutectCalls:
#     input:
#         vcf=output_somatic+"{sample}_somatic_mutect2.vcf",
#         idx=output_somatic+"{sample}_somatic_mutect2.vcf.idx",
#         stats=output_somatic+"{sample}_somatic_mutect2.vcf.stats",
#         vcf_stats=output_somatic+"{sample}_somatic_mutect2.vcf.stats",
#         contamination=output_somatic+"{sample}_merged_contamination.table",
#         #read_orientation=output_somatic+"{sample}_read-orientation-model.tar.gz"
#     output:
#         vcf=output_somatic+"{sample}_somatic_mutect2_filtered.vcf",
#         idx=output_somatic+"{sample}_somatic_mutect2_filtered.vcf.idx"
#     shell:
#         ""
#         gatk --java-options {mem} FilterMutectCalls \
#         -R {ref} \
#         -V {input.vcf} \
#         --contamination-table {input.contamination} \
#         --stats {input.stats} \
#         -L {interval_list} \
#         -O {output.vcf}
#         ""
#         # --orientation-bias-artifact-priors {input.read_orientation} \

rule FilterMutectCalls:
    input:
        vcf=output_somatic+"{sample}_somatic_mutect2.vcf",
        idx=output_somatic+"{sample}_somatic_mutect2.vcf.idx",
        stats=output_somatic+"{sample}_somatic_mutect2.vcf.stats",
        vcf_stats=output_somatic+"{sample}_somatic_mutect2.vcf.stats",
        contamination=output_somatic+"{sample}_merged_contamination.table",
        read_orientation=output_somatic+"{sample}_read-orientation-model.tar.gz"
    output:
        vcf=output_somatic+"{sample}_somatic_mutect2_filtered.vcf",
        idx=output_somatic+"{sample}_somatic_mutect2_filtered.vcf.idx"
    shell:
        """
        gatk --java-options {mem} FilterMutectCalls \
        -R {ref} \
        -V {input.vcf} \
        --contamination-table {input.contamination} \
        --orientation-bias-artifact-priors {input.read_orientation} \
        --stats {input.stats} \
        -L {interval_list} \
        -O {output.vcf}
        """
