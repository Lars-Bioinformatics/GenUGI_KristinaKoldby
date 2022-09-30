__title__ = "Pipeline for Somatic JOINT Variant Calling with Varscan2 - Kristina's project"
__author__ = "Lars Andersen <larsmew@gmail.com> and Kristina Koldby"
__date__ = "8/1/2020"
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
target_regions = "/work/sdukoldby/resources/hg38/target_regions/bed_files/"

# SAMPLES = ["ECV2-4"]
# TUMORS = ["ECV2-4-biopsi-H1_tumor","ECV2-4-biopsi-H2_tumor","ECV2-4-op-1401_tumor","ECV2-4-op-1801_tumor","ECV2-4-op-1802_tumor","ECV2-4-plasma170503_tumor","ECV2-4-plasma170712_tumor","ECV2-4-plasma180205_tumor"]
# NORMALS = ["ECV2-4-blod_normal"]
# SAMPLES = ["ECV2-8"]
# TUMORS = ["ECV2-8-biopsi-I1_tumor","ECV2-8-biopsi-I2_tumor","ECV2-8-op-01_tumor","ECV2-8-op-02_tumor","ECV2-8-op-03_tumor","ECV2-8-plasma170522_tumor","ECV2-8-plasma170728_tumor","ECV2-8-plasma180806_tumor"]
# NORMALS = ["ECV2-8-blod_normal"]
# SAMPLES = ["ECV2-29"]
# TUMORS = ["ECV2-29-biopsi-C1_tumor","ECV2-29-biopsi-C2_tumor","ECV2-29-op-A1_tumor","ECV2-29-op-A2_tumor","ECV2-29-op-B1_tumor","ECV2-29-op-B2_tumor","ECV2-29-plasma171124_tumor","ECV2-29-plasma180119_tumor","ECV2-29-plasma180619_tumor"]
# NORMALS = ["ECV2-29-blod_normal"]
# SAMPLES = ["ECV2-31"]
# TUMORS = ["ECV2-31-biopsi-F1_tumor","ECV2-31-biopsi-F2_tumor","ECV2-31-op-D1_tumor","ECV2-31-op-D2_tumor","ECV2-31-op-E1_tumor","ECV2-31-op-E2_tumor","ECV2-31-plasma171215_tumor","ECV2-31-plasma180202_tumor"]
# NORMALS = ["ECV2-31-blod_normal"]
# SAMPLES = ["ECV2-35"]
# TUMORS = ["ECV2-35-biopsi-G1_tumor","ECV2-35-biopsi-G2_tumor","ECV2-35-op-01_tumor","ECV2-35-op-02_tumor","ECV2-35-op-03_tumor","ECV2-35-op-04_tumor","ECV2-35-op-06_tumor","ECV2-35-op-07_tumor","ECV2-35-plasma180102_tumor","ECV2-35-plasma180316_tumor","ECV2-35-plasma180601_tumor"]
# NORMALS = ["ECV2-35-blod_normal"]
# SAMPLES = ["PC1-10"]
# TUMORS = ["PC1-10-op-P1_tumor","PC1-10-op-P2_tumor","PC1-10-op-P3_tumor","PC1-10-op-P4_tumor","PC1-10-plasma170928_tumor"]
# NORMALS = ["PC1-10-blod_normal"]
# SAMPLES = ["PC1-14"]
# TUMORS = ["PC1-14-op-P5_tumor","PC1-14-op-P6_tumor","PC1-14-op-P7_tumor","PC1-14-op-P8_tumor","PC1-14-plasma180111_tumor","PC1-14-plasma180907_tumor"]
# NORMALS = ["PC1-14-blod_normal"]
# SAMPLES = ["PC1-18"]
# TUMORS = ["PC1-18-op-P9_tumor","PC1-18-op-P10_tumor","PC1-18-op-P11_tumor","PC1-18-op-P12_tumor","PC1-18-plasma180219_tumor"]
# NORMALS = ["PC1-18-blod_normal"]

# print(SAMPLES)
# print(TUMORS)
# print(NORMALS)

PATIENTS = ["ECV2-4", "ECV2-8", "ECV2-29", "ECV2-31", "ECV2-35", "PC1-10", "PC1-14", "PC1-18", "PON-all"]
# PATIENTS = ["ECV2-4"]


INTERVALS, = glob_wildcards(target_regions+"{interval}.bed")
INTERVALS = sorted(INTERVALS)
print(INTERVALS)
print (PATIENTS)


#########################################################
####                      Output                     ####
#########################################################
output_somatic = "varscan_somatic_joint_calling/"


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
    shell("mkdir -p "+output_somatic)
    shell("mkdir -p "+output_somatic+"split")
    shell("mkdir -p "+output_somatic+"split_vcf")


#########################################################
####                  Run All Rules                  ####
#########################################################
'''
Rule all
'''
print(expand(output_somatic+"{sample}_cns_varscan2.vcf.gz", sample=PATIENTS))
# sys.exit()

rule all_pairs:
    input:
        expand(output_somatic+"{sample}_cns_varscan2.vcf.gz", sample=PATIENTS),
        expand(output_somatic+"{sample}_cns_varscan2_somatic_pon-filtered.vcf.gz", sample=PATIENTS),
        expand(output_somatic+"{sample}_cns_varscan2_somatic_pon-filtered_strict.vcf.gz", sample=PATIENTS),
        # expand(output_somatic+"{sample}_cns_varscan2_somatic.vcf.gz", sample=PATIENTS),
        # expand(output_somatic+"{sample}_cns_varscan2_somatic_NormADzero.vcf.gz", sample=PATIENTS),
        # expand(output_somatic+"{sample}_sampleList.txt", sample=PATIENTS)


##########################################################
####  Create list of sample names                     ####
##########################################################
rule sampleList:
    output:
        output_somatic+"{sample}_sampleList.txt"
    params:
        sample_names = lambda wildcards: config[wildcards.sample]["all"]
    run:
        # print("_tagseq-medexome\n".join(SAMPLE_NAMES))
        # print(output[0])
        with open(output[0],"w") as f:
            f.write("\n".join([s for s in params[0]]))
            # f.write("\n".join([s+"_tagseq-medexome" for s in params[0]]))


##########################################################
####  Samtools mpileup                                ####
##########################################################
rule mpileup:
    input:
        bam = lambda wildcards: [b+".connor.recalibrated.bam" for b in config[wildcards.sample]["all"]],
        intervals=target_regions+"{interval}.bed"
    output:
        mpileups=output_somatic+"split/{sample}__{interval}__split.mpileup"
    shell:
        """
        samtools mpileup \
        -f {ref} \
        -B \
        -l {input.intervals}\
        -d 1000000 \
        {input.bam} > {output}
        """


##########################################################
####  Call Somatic Variants using Varscan2 on matched ####
####   Tumor-Normal samples                           ####
##########################################################
'''
Multi-sample calling using Varscan2
'''
rule Varscan2:
    input:
        mpileups=output_somatic+"split/{sample}__{interval}__split.mpileup",
        sample_list=output_somatic+"{sample}_sampleList.txt"
    output:
        vcf=output_somatic+"split_vcf/{sample}_cns__{interval}__split_varscan2.vcf"
    shell:
        """
        varscan mpileup2cns {input.mpileups} \
        --min-coverage 1 \
        --min-reads2 2 \
        --min-var-freq 0 \
        --output-vcf 1 \
        --min-avg-qual 20 \
        --variants 1 \
        --strand-filter 0 \
        --vcf-sample-list {input.sample_list} \
        > {output}
        """

# rule index_vcf:
#     input:
#         vcf_subfile=output_somatic+"split_vcf/{sample}_cns__{interval}__split_varscan2.vcf"
#     output:
#         vcf_index=output_somatic+"split_vcf/{sample}_cns__{interval}__split_varscan2.vcf.idx"
#     shell:
#         """
#         gatk IndexFeatureFile -I {input}
#         """
#
#
# rule merge_Varscan2:
#     input:
#         vcf_subfile=expand(output_somatic+"split_vcf/{{sample}}_cns__{interval}__split_varscan2.vcf", interval=INTERVALS)
#     output:
#         vcf=output_somatic+"{sample}_cns_varscan2.vcf"
#     params:
#         vcf_subfile=expand("-I "+output_somatic+"split_vcf/{{sample}}_cns__{interval}__split_varscan2.vcf", interval=INTERVALS)
#     shell:
#         """
#         gatk --java-options {mem} GatherVcfs \
#         {params.vcf_subfile} \
#         -O {output.vcf}
#         """

rule index_vcf:
    input:
        vcf_subfile=output_somatic+"split_vcf/{sample}_cns__{interval}__split_varscan2.vcf"
    output:
        vcf_subfile_compressed=output_somatic+"split_vcf/{sample}_cns__{interval}__split_varscan2.vcf.gz",
        vcf_index=output_somatic+"split_vcf/{sample}_cns__{interval}__split_varscan2.vcf.gz.tbi"
    shell:
        """
        bgzip -c {input} > {output.vcf_subfile_compressed}
        bcftools index {output.vcf_subfile_compressed} -o {output.vcf_index}
        """

rule merge_Varscan2:
    input:
        vcf_subfile=expand(output_somatic+"split_vcf/{{sample}}_cns__{interval}__split_varscan2.vcf.gz", interval=INTERVALS),
        vcf_index=expand(output_somatic+"split_vcf/{{sample}}_cns__{interval}__split_varscan2.vcf.gz.tbi", interval=INTERVALS)
    output:
        vcf=output_somatic+"{sample}_cns_varscan2.vcf.gz",
        vcf_index=output_somatic+"{sample}_cns_varscan2.vcf.gz.tbi"
    params:
        temp=output_somatic+"{sample}_cns_varscan2.vcf"
    shell:
        """
        bcftools concat \
        {input.vcf_subfile} \
        --no-version \
        -o {params.temp};
        bgzip {params.temp};
        tabix -p vcf {output.vcf}
        """

        # """
        # bcftools concat \
        # {input.vcf_subfile} \
        # --no-version \
        # -o {output.vcf};
        # bgzip -c {output.vcf} > {output.vcf}.gz;
        # tabix -p vcf {output.vcf}.gz
        # """

### Remove variants found with more than one read i blood sample ###
### IMPORTANT: 'FORMAT/AD[0]' requires that blood sample is the first sample in joint vcf
### This is the case when blood sample is listed as first sample in category "all" i config file

rule somatic_filter:
    input:
        output_somatic+"{sample}_cns_varscan2.vcf.gz"
    output:
        output_somatic+"{sample}_cns_varscan2_somatic.vcf.gz"
    shell:
        """
        bcftools filter \
        --include 'FORMAT/AD[0] <= 1' \
        {input} \
        -o {output}
        """

rule somatic_filter_strict:
    input:
        output_somatic+"{sample}_cns_varscan2.vcf.gz"
    output:
        output_somatic+"{sample}_cns_varscan2_somatic_NormADzero.vcf.gz"
    shell:
        """
        bcftools filter \
        --include 'FORMAT/AD[0] = 0' \
        {input} \
        -o {output}
        """

rule pon_filtering:
    input:
        vcf = output_somatic+"{sample}_cns_varscan2.vcf.gz",
        pon = output_somatic+"PON-all_cns_varscan2.vcf.gz"
    output:
        output_somatic+"{sample}_cns_varscan2_somatic_pon-filtered.vcf.gz"
    shell:
        """
        bcftools isec \
        -C \
        -c none \
        -w1 \
        {input.vcf} \
        {input.pon} \
        -o {output}
        """

rule pon_filtering_strict:
    input:
        vcf = output_somatic+"{sample}_cns_varscan2.vcf.gz",
        pon = output_somatic+"PON-all_cns_varscan2.vcf.gz"
    output:
        output_somatic+"{sample}_cns_varscan2_somatic_pon-filtered_strict.vcf.gz"
    shell:
        """
        bcftools isec \
        -C \
        -c none \
        -w1 \
        {input.vcf} \
        {input.pon} | \
        bcftools filter \
        --include 'FORMAT/AD[0] = 0' \
        -o {output}
        """
