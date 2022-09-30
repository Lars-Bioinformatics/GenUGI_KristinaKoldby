## Run from folder with vcf files
## No need to rearrange sample order in vcf files if created using "split_joint_vcfs.sh" - in that case tumor is always first and blood second

#### REMEMBER TO ADJUST rule filter_vcf ACCORDING TO VARIANT-CALLER
#### REMEMBER TO ADJUST rule tmp_file1 ACCORDING TO VARIANT-CALLER

window_size = 1000
# variant_caller = "mutect2"
variant_caller = "varscan2"
ADmin = 20
DPmin = 20

# bam_path = "/work/sdukoldby/data/G45-2016_genugi/exome_fastq/exome_fastq_merged_2nd_novaseq_run/hg38/Connor/"
# ngcgh_path = "/work/sdukoldby/data/G45-2016_genugi/exome_fastq/exome_fastq_merged_2nd_novaseq_run/hg38/Connor/ngcgh/"
# ngcgh_sorted_path = "/work/sdukoldby/data/G45-2016_genugi/exome_fastq/exome_fastq_merged_2nd_novaseq_run/hg38/Connor/ngcgh_sorted/"
bam_path = "/work/sdukoldby/data/G45-2016_genugi/exome_fastq/exome_fastq_merged_2nd_novaseq_run/hg38/sample_merge/Connor/"
ngcgh_path = "/work/sdukoldby/data/G45-2016_genugi/exome_fastq/exome_fastq_merged_2nd_novaseq_run/hg38/sample_merge/Connor/ngcgh/"
ngcgh_sorted_path = "/work/sdukoldby/data/G45-2016_genugi/exome_fastq/exome_fastq_merged_2nd_novaseq_run/hg38/sample_merge/Connor/ngcgh_sorted/"
header = "vcf_header.txt"

SAMPLES, = glob_wildcards("{sample}.vcf.gz")
# SAMPLES = "ECV2-4-biopsi-H1_tumor_tagseq-medexome-deep-seq"
# SAMPLES, = glob_wildcards("{sample}_tumor_tagseq-medexome-deep-seq.vcf.gz","{sample}_tumor_tagseq-medexome.vcf.gz")
TUMORS, = glob_wildcards(ngcgh_path + "{tumor}_medexome-deep-seq_ngcgh_w1000.txt")
# TUMORS = "ECV2-4-biopsi-H1_tumor"
#SAMPLES = "ECV2-4-biopsi-H1_tumor_vs_ECV2-4-blod"
# PAIR, = glob_wildcards("{sample}.vcf.gz")
# PAIR = ["ECV2-29-biopsi-C1_tumor_tagseq-medexome-deep-seq"]
# TUMORS, = lambda wildcards: config[wildcards.sample]["tumor"]]

# configfile: "samples_tumor.yaml"
configfile: "samples_vcf.yaml"

# print(SAMPLES)
# print(TUMORS)
# sys.exit()

onstart:
    shell("mkdir -p vaf_log2r_plots/"),
    shell("mkdir -p vaf_log2r_plots/plot_tables/"),
    shell("mkdir -p vaf_log2r_plots/plot_tables/ADmin{ADmin}")

rule all:
    input:
        expand("vaf_log2r_plots/plot_tables/ADmin{ADmin}/{sample}_somatic_{variant_caller}_log2r_w{window_size}_af_filtered_ADmin{ADmin}-DPmin{DPmin}.txt",
            sample=SAMPLES,
            variant_caller=variant_caller,
            ADmin=ADmin,
            DPmin=DPmin,
            window_size=window_size)


#################################################################
############# Add Log2R from ngcgh file to vcf file #############
#################################################################

#### Sort, bgzip and index ngcgh-file (required by bcftools) ####

rule sort_ngcgh:
    input:
        ngcgh = ngcgh_path + "{tumor}_medexome-deep-seq_ngcgh_w{window_size}.txt"
    output:
        temp(ngcgh_sorted_path + "{tumor}_medexome-deep-seq_ngcgh_w{window_size}_sorted.txt")
    shell:
        """
        sort -k1,1 -k2,2n -k3,3n {input} > {output}
        """
# ngcgh_path + "{tumor}_medexome-deep-seq_ngcgh_w{window_size}.txt"
# lambda wildcards: [s+"_medexome-deep-seq_ngcgh_w{window_size}.txt" for s in config[wildcards.sample]["tumor"]]

rule make_header:
    input:
        ngcgh_sorted_path + "{tumor}_medexome-deep-seq_ngcgh_w{window_size}_sorted.txt"
    output:
        temp(ngcgh_sorted_path + "{tumor}_medexome-deep-seq_ngcgh_w{window_size}_sorted_header.txt")
    shell:
        """
        echo -e "CHROM\tFROM\tTO\tNORMCOUNT\tTUMORCOUNT\tLog2R" | \
        cat - {input} > {output}
        """

rule bgzip_ngcgh:
    input:
        ngcgh_sorted_path + "{tumor}_medexome-deep-seq_ngcgh_w{window_size}_sorted_header.txt"
    output:
        ngcgh_sorted_path + "{tumor}_medexome-deep-seq_ngcgh_w{window_size}_sorted_header.txt.gz"
    shell:
        """
        bgzip {input}
        """

rule index_ngcgh:
    input:
        ngcgh_sorted_path + "{tumor}_medexome-deep-seq_ngcgh_w{window_size}_sorted_header.txt.gz"
    output:
        ngcgh_sorted_path + "{tumor}_medexome-deep-seq_ngcgh_w{window_size}_sorted_header.txt.gz.tbi"
    shell:
        """
        tabix -s1 -b2 -e3 -S1 {input}
        """


#### Append log2-ratio from ngCGH-file ####

rule add_log2r_to_vcf_1:
    input:
        ngcgh = ngcgh_sorted_path + "{tumor}_medexome-deep-seq_ngcgh_w{window_size}_sorted_header.txt.gz",
        index = ngcgh_sorted_path + "{tumor}_medexome-deep-seq_ngcgh_w{window_size}_sorted_header.txt.gz.tbi",
        vcf = "{tumor}_tagseq-medexome-deep-seq.vcf.gz",
    output:
        "vaf_log2r_plots/{tumor}_tagseq-medexome-deep-seq_somatic_{variant_caller}_log2r_w{window_size}.vcf"
    shell:
        """
        bcftools annotate \
        -a {input.ngcgh} \
        -c CHROM,FROM,TO,-,-,LOG2R \
        -h {header} \
        {input.vcf} \
        -o {output}
        """

rule add_log2r_to_vcf_2:
    input:
        ngcgh = ngcgh_sorted_path + "{tumor}_medexome-deep-seq_ngcgh_w{window_size}_sorted_header.txt.gz",
        index = ngcgh_sorted_path + "{tumor}_medexome-deep-seq_ngcgh_w{window_size}_sorted_header.txt.gz.tbi",
        vcf = "{tumor}_tagseq-medexome.vcf.gz",
    output:
        "vaf_log2r_plots/{tumor}_tagseq-medexome_somatic_{variant_caller}_log2r_w{window_size}.vcf"
    shell:
        """
        bcftools annotate \
        -a {input.ngcgh} \
        -c CHROM,FROM,TO,-,-,LOG2R \
        -h {header} \
        {input.vcf} \
        -o {output}
        """

# input:
        # ngcgh = lambda wildcards: ["vaf_log2r_plots/ngcgh_sorted/"+s+"_medexome-deep-seq_ngcgh_w{window_size}_sorted_header.txt.gz" for s in config[wildcards.sample]['tumor']],
        # index = lambda wildcards: ["vaf_log2r_plots/ngcgh_sorted/"+s+"_medexome-deep-seq_ngcgh_w{window_size}_sorted_header.txt.gz.tbi" for s in config[wildcards.sample]['tumor']],
        # vcf = "{sample}.vcf.gz",
        # ngcgh = "vaf_log2r_plots/ngcgh_sorted/"+[b+"_tumor_medexome-deep-seq_ngcgh_w{window_size}_sorted_header.txt.gz" for b in config["tumor"][wildcards.sample]],
        # index = "vaf_log2r_plots/ngcgh_sorted/"+[b+"_tumor_medexome-deep-seq_ngcgh_w{window_size}_sorted_header.txt.gz.tbi" for b in config["tumor"][wildcards.sample]],
        # vcf = "{sample}.vcf.gz",
        # ngcgh = "vaf_log2r_plots/ngcgh_sorted/{tumor}_medexome-deep-seq_ngcgh_w{window_size}_sorted_header.txt.gz",
        # index = "vaf_log2r_plots/ngcgh_sorted/{tumor}_medexome-deep-seq_ngcgh_w{window_size}_sorted_header.txt.gz.tbi",
        # vcf = lambda wildcards: [b+".vcf.gz" for b in config["vcf"][wildcards.sample]],
        # ngcgh = "vaf_log2r_plots/ngcgh_sorted/{tumor}_medexome-deep-seq_ngcgh_w{window_size}_sorted_header.txt.gz", tumor=config[tumor],
        # index = "vaf_log2r_plots/ngcgh_sorted/{tumor}_medexome-deep-seq_ngcgh_w{window_size}_sorted_header.txt.gz.tbi",
        # vcf = "{sample}.vcf.gz",

########################################################################
######### LAV FIL TIL PLOTS (AF AND LOG2R ON "PASS" VARIANTS) ##########
########################################################################

#### Filter vcf file on number of variant reads ####
#### "FORMAT/AD[0:1]" refers to the AD value in the format column, first sample, second value, ####
#### which corresponds to the AD value of the variant allele in the tumor sample ####
#### "FORMAT/DP[1]" refers to the DP value in the format column, second sample, i.e. depth in normal sample ####


rule filter_vcf:
    input:
        "vaf_log2r_plots/{sample}_somatic_{variant_caller}_log2r_w{window_size}.vcf"
    output:
        "vaf_log2r_plots/{sample}_somatic_{variant_caller}_log2r_w{window_size}_filtered_ADmin{ADmin}-DPmin{DPmin}.vcf"
    shell:
        """
        bcftools filter \
        --include 'FORMAT/AD[0] >= {ADmin} && FORMAT/DP[1] >= {DPmin}' \
        {input} \
        -o {output}
        """

### MUTECT USE THIS ###
### --include 'FORMAT/AD[0:1] >= {ADmin} && FORMAT/DP[1] >= {DPmin}' \

### VARSCAN USE THIS ###
### --include 'FORMAT/AD[0] >= {ADmin} && FORMAT/DP[1] >= {DPmin}' \


#### Create new file with values from columns #CHROM ($1), POS ($2), LOG2R ($8) and AF in tumor-sample ($11) ####

rule tmp_file1:
    input:
        "vaf_log2r_plots/{sample}_somatic_{variant_caller}_log2r_w{window_size}_filtered_ADmin{ADmin}-DPmin{DPmin}.vcf"
    output:
        temp("vaf_log2r_plots/{sample}_somatic_{variant_caller}_log2r_w{window_size}_passed_ADmin{ADmin}-DPmin{DPmin}_tmp1.txt")
    shell:
        """
        bcftools query \
        -H \
        -s {wildcards.sample} \
        --include 'FILTER = "PASS"' \
        -f '%CHROM\t%POS\t[%FREQ]\t%INFO/LOG2R\n' \
        {input} \
        -o {output}
        """

### MUTECT USE THIS ###
### -f '%CHROM\t%POS\t[%AF]\t%INFO/LOG2R\n' \

### VARSCAN USE THIS ###
### -f '%CHROM\t%POS\t[%FREQ]\t%INFO/LOG2R\n' \


#### Remove "LOG2R=" from file

rule tmp_file2:
    input:
        "vaf_log2r_plots/{sample}_somatic_{variant_caller}_log2r_w{window_size}_passed_ADmin{ADmin}-DPmin{DPmin}_tmp1.txt"
    output:
        temp("vaf_log2r_plots/{sample}_somatic_{variant_caller}_log2r_w{window_size}_passed_ADmin{ADmin}-DPmin{DPmin}_tmp2.txt")
    shell:
        """
        sed 's/LOG2R=//g' {input} > {output}
        """

#### Add headers to file

rule final_plot_file:
    input:
        "vaf_log2r_plots/{sample}_somatic_{variant_caller}_log2r_w{window_size}_passed_ADmin{ADmin}-DPmin{DPmin}_tmp2.txt"
    output:
        "vaf_log2r_plots/plot_tables/ADmin{ADmin}/{sample}_somatic_{variant_caller}_log2r_w{window_size}_af_filtered_ADmin{ADmin}-DPmin{DPmin}.txt"
    shell:
        """
        echo -e "CHR\tPOS\tAF\tLOG2R" | \
        cat - {input} | \
        awk -v OFS="\t" '$1=$1' - \
        > {output}
        """

###############################
##### Create plots (in R) #####
###############################

# module load R
# R
#
# files <- list.files(path=".", pattern="*.txt", full.names=TRUE, recursive=FALSE)
# lapply(files, function(x) {
#
# data<-as.matrix(read.table(x,header=TRUE, sep="\t"))
#
# xval<-as.numeric(sub("%","",data[,3]))/100
# yval<-as.numeric(data[,4])
#
# png(paste(strsplit(x,".txt")[[1]],".png",sep=""))
#
# plot(xval,yval,xlab="VAF",ylab="Log2R",xlim=c(0,1),ylim=c(-3,3),pch=c(rep(18,length(xval))))
# abline(h=0,lwd=2)
#
# dev.off()
#
# rm(list=ls())
#
# })
#
# ### Mutect use:
# xval<-as.numeric(data[,3])
# yval<-as.numeric(data[,4])
#
# ### Varscan use:
# xval<-as.numeric(sub("%","",data[,3]))/100
# yval<-as.numeric(data[,4])
