## Run from folder with vcf files ("paired_tumor_samples_vcfs")
## No need to rearrange sample order in vcf files if created using "split_joint_vcfs.sh" - in that case tumor is always first and blood second

###################################################################
#### REMEMBER TO ADJUST BOTH RULES ACCORDING TO VARIANT-CALLER ####
###################################################################

# window_size = 1000
# variant_caller = "mutect2"
# variant_caller = "varscan2"
ADmin = 10
DPmin = 20

SAMPLES, = glob_wildcards("{sample}.vcf.gz")
# SAMPLES = "ECV2-29-biopsi-C1_vs_ECV2-29-biopsi-C2_somatic_mutect2"

onstart:
    shell("mkdir -p vaf_vaf_plots/"),
    shell("mkdir -p vaf_vaf_plots/ADmin{ADmin}"),
    # shell("mkdir -p vaf_vaf_plots/plot_tables/")


rule all:
    input:
        expand("vaf_vaf_plots/ADmin{ADmin}/{sample}_filtered_ADmin{ADmin}-DPmin{DPmin}-pass_plot-table.txt",
            sample=SAMPLES,
            ADmin=ADmin,
            DPmin=DPmin)

# variant_caller=variant_caller,
# window_size=window_size,

##################################################################
######### LAV FIL TIL PLOTS: VAF-sample1 vs VAF-sample2 ##########
##################################################################

#### Filter vcf file on number of variant reads ####
#### "FORMAT/AD[*:1]" refers to the AD value in the format column, in any of the one of the samples, second value (alt allel), ####
#### Requires that vcf has been filtered beforehand and only contains somatic variants ####
#### "FORMAT/DP[2]" refers to the DP value in the format column, third sample, i.e. depth in normal sample ####


rule filter_vcf:
    input:
        "{sample}.vcf.gz"
    output:
        "vaf_vaf_plots/ADmin{ADmin}/{sample}_filtered_ADmin{ADmin}-DPmin{DPmin}-pass.vcf"
    shell:
        """
        bcftools filter \
        --include 'FORMAT/AD[*:1] >= {ADmin} && FORMAT/DP[2] >= {DPmin} && FILTER = "PASS"' \
        {input} \
        -o {output}
        """

### MUTECT USE THIS ###
### --include 'FORMAT/AD[*:1] >= {ADmin} && FORMAT/DP[2] >= {DPmin}' \

### VARSCAN USE THIS ###
### --include 'FORMAT/AD[*] >= {ADmin} && FORMAT/DP[2] >= {DPmin}' \


rule plot_table:
    input:
        "vaf_vaf_plots/ADmin{ADmin}/{sample}_filtered_ADmin{ADmin}-DPmin{DPmin}-pass.vcf"
    output:
        "vaf_vaf_plots/ADmin{ADmin}/{sample}_filtered_ADmin{ADmin}-DPmin{DPmin}-pass_plot-table.txt"
    shell:
        """
        bcftools query \
        -H \
        -f '%CHROM\t%POS\t[%AF\t]\n' \
        {input} | \
        sed 's/LOG2R=//g' > \
        {output}
        """

### MUTECT USE THIS ###
### -f '%CHROM\t%POS\t[%AF\t]\n' \

### VARSCAN USE THIS ###
### -f '%CHROM\t%POS\t[%FREQ\t]\n' \

# -s {wildcards.sample} \
# --include 'FILTER = "PASS"' \


###############################
##### Create plots (in R) #####
###############################

# module load R
# R
#
# files <- list.files(path=".", pattern="*_plot-table.txt", full.names=TRUE, recursive=FALSE)
# lapply(files, function(x) {
#
# data<-as.matrix(read.table(x,header=TRUE, sep="\t"))
#
# xval<-as.numeric(data[,4])
# yval<-as.numeric(data[,3])
#
# png(paste(strsplit(x,"_plot-table.txt")[[1]],".png",sep=""))
#
# plot(xval,yval,xlab="VAF sample 2",ylab="VAF sample 1",xlim=c(0,1),ylim=c(0,1),pch=c(rep(18,length(xval))))
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
# yval<-as.numeric(sub("%","",data[,4]))/100
