## Run from folder with vcf files
######## RUN THESE COMMANDS BEFORE SCRIPT ########
## mkdir -p original_biopsy_vcf
## mv *biopsi*_filtered.vcf original_biopsy_vcf/
## cd original_biopsy_vcf/
## rename mutect2_filtered mutect2_original_filtered *
## bash rearrange_vcfs.sh

SAMPLES, = glob_wildcards("{sample}_normal_tagseq-medexome_somatic_mutect2_filtered.vcf")
#SAMPLES = "ECV2-4-biopsi-H1_tumor_vs_ECV2-4-blod"

header = "vcf_header.txt"
ngcgh_path = "/work/sdukoldby/data/G45-2016_genugi/exome_fastq/exome_fastq_merged/hg38/Connor/ngcgh_lars/"

onstart:
    shell("mkdir -p vaf_log2r_plots/"),
    shell("mkdir -p vaf_log2r_plots/plot_tables/")

rule all:
    input:
        expand("vaf_log2r_plots/plot_tables/{sample}_somatic_mutect2_log2r_af_filtered_ADmin4-DPmin50.txt", sample=SAMPLES)


#################################################################
############# Add Log2R from ngcgh file to vcf file #############
#################################################################

#### Sort, bgzip and index ngcgh-file (required by bcftools) ####

rule sort_ngcgh:
    input:
        ngcgh_path + "{sample}_normal_tagseq-medexome_ngcgh_w1000.txt"
    output:
        temp("vaf_log2r_plots/{sample}_normal_tagseq-medexome_ngcgh_w1000_sorted.txt")
    shell:
        """
        sort -k1,1 -k2,2n -k3,3n {input} > {output}
        """

rule make_header:
    input:
        "vaf_log2r_plots/{sample}_normal_tagseq-medexome_ngcgh_w1000_sorted.txt"
    output:
        temp("vaf_log2r_plots/{sample}_normal_tagseq-medexome_ngcgh_w1000_sorted_header.txt")
    shell:
        """
        echo -e "CHROM\tFROM\tTO\tNORMCOUNT\tTUMORCOUNT\tLog2R" | \
        cat - {input} > {output}
        """

rule bgzip_ngcgh:
    input:
        "vaf_log2r_plots/{sample}_normal_tagseq-medexome_ngcgh_w1000_sorted_header.txt"
    output:
        "vaf_log2r_plots/{sample}_normal_tagseq-medexome_ngcgh_w1000_sorted_header.txt.gz"
    shell:
        """
        bgzip {input}
        """

rule index_ngcgh:
    input:
        "vaf_log2r_plots/{sample}_normal_tagseq-medexome_ngcgh_w1000_sorted_header.txt.gz"
    output:
        "vaf_log2r_plots/{sample}_normal_tagseq-medexome_ngcgh_w1000_sorted_header.txt.gz.tbi"
    shell:
        """
        tabix -s1 -b2 -e3 -S1 {input}
        """


#### Append log2-ratio from ngCGH-file ####

rule add_log2r_to_vcf:
    input:
        ngcgh = "vaf_log2r_plots/{sample}_normal_tagseq-medexome_ngcgh_w1000_sorted_header.txt.gz",
        index = "vaf_log2r_plots/{sample}_normal_tagseq-medexome_ngcgh_w1000_sorted_header.txt.gz.tbi",
        vcf = "{sample}_normal_tagseq-medexome_somatic_mutect2_filtered.vcf"
    output:
        "vaf_log2r_plots/{sample}_somatic_mutect2_log2r.vcf"
    shell:
        """
        bcftools annotate \
        -a {input.ngcgh} \
        -c CHROM,FROM,TO,-,-,LOG2R \
        -h {header} \
        {input.vcf} \
        -o {output}
        """


########################################################################
######### LAV FIL TIL PLOTS (AF AND LOG2R ON "PASS" VARIANTS) ##########
########################################################################

#### Filter vcf file on number of variant reads ####
#### "FORMAT/AD[1:1]" refers to the AD value in the format column, second sample, second value, ####
#### which corresponds to the AD value of the variant allele in the tumor sample ####

rule filter_vcf:
    input:
        "vaf_log2r_plots/{sample}_somatic_mutect2_log2r.vcf"
    output:
        "vaf_log2r_plots/{sample}_somatic_mutect2_log2r_filtered_ADmin4-DPmin50.vcf"
    shell:
        """
        bcftools filter \
        --include 'FORMAT/AD[1:1] > 3 && FORMAT/DP[0] >= 50' \
        {input} \
        -o {output}
        """

## && FORMAT/DP[0] >= 20

#### Create file with only variants that has FILTER=PASS ####
#### Behold overskrifter på kolonner ####

rule passed_variants:
    input:
        "vaf_log2r_plots/{sample}_somatic_mutect2_log2r_filtered_ADmin4-DPmin50.vcf"
    output:
        "vaf_log2r_plots/{sample}_somatic_mutect2_log2r_filtered_ADmin4-DPmin50_passed.txt"
    shell:
        """
        awk \
        -v FS='\t' \
        -v OFS='\t' \
        '{{if($7=="PASS" || $7=="FILTER")print $0}}' \
        {input} > {output}
        """


#### Create new file with values from columns #CHROM ($1), POS ($2), LOG2R ($8) and AF in tumor-sample ($11) ####

rule tmp_file1:
    input:
        "vaf_log2r_plots/{sample}_somatic_mutect2_log2r_filtered_ADmin4-DPmin50_passed.txt"
    output:
        temp("vaf_log2r_plots/{sample}_somatic_mutect2_log2r_passed_tmp1.txt")
    shell:
        """
        awk -v FS='\t' -v OFS=':' \
        '{{print $1 ":" $2 ":" $11 ":" $8}}' \
        {input} | \
        awk -v FS=':' -v OFS='\t' \
        '{{$3=$4=""; print $0}}' - | \
        awk  '{{gsub(";","\t",$0); print $1, $2, $3, $NF;}}' - | \
        grep 'LOG2R' > {output}
        """

#### Remove "LOG2R=" from file

rule tmp_file2:
    input:
        "vaf_log2r_plots/{sample}_somatic_mutect2_log2r_passed_tmp1.txt"
    output:
        temp("vaf_log2r_plots/{sample}_somatic_mutect2_log2r_passed_tmp2.txt")
    shell:
        """
        sed 's/LOG2R=//g' {input} > {output}
        """

#### Add headers to file

rule final_plot_file:
    input:
        "vaf_log2r_plots/{sample}_somatic_mutect2_log2r_passed_tmp2.txt"
    output:
        "vaf_log2r_plots/plot_tables/{sample}_somatic_mutect2_log2r_af_filtered_ADmin4-DPmin50.txt"
    shell:
        """
        echo -e "CHR\tPOS\tAF\tLOG2R" | \
        cat - {input} | \
        awk -v OFS="\t" '$1=$1' - \
        > {output}
        """
