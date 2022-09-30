## Run from folder with vcf files

SAMPLES, = glob_wildcards("{sample}_cns_varscan2.vcf.gz")
# SAMPLES, = "ECV2-4"

configfile: "/work/sdukoldby/data/G45-2016_genugi/exome_fastq/exome_fastq_merged/hg38/Connor/samples.yaml"

header = "vcf_header.txt"
ngcgh_path = "/work/sdukoldby/data/G45-2016_genugi/exome_fastq/exome_fastq_merged/hg38/Connor/ngcgh_lars/"

onstart:
    shell("mkdir -p vaf_log2r_plots/"),
    shell("mkdir -p vaf_log2r_plots/plot_tables/")
    shell("mkdir -p vaf_log2r_plots/private_variants/")
    shell("mkdir -p vaf_log2r_plots/private_variants/plot_tables/")

rule all:
    input:
        # expand("vaf_log2r_plots/{sample}_varscan_filtered-NormADmax1-TumorADmin7-NormDPmin20.vcf.gz", sample=SAMPLES)
        #expand("vaf_log2r_plots/plot_tables/{sample}_somatic_varscan_log2r_af_filtered_ADmin3.txt", sample=SAMPLES)
        # [expand("vaf_log2r_plots/{sample_subfile}_single_varscan2.vcf.gz", sample_subfile=config[sample]["tumors"]) for sample in SAMPLES]
        # [expand("vaf_log2r_plots/{sample_subfile}_single_varscan_filtered-NormADmax1-TumorADmin7-NormDPmin20.vcf", sample_subfile=config[sample]["tumors"]) for sample in SAMPLES]
        # [expand("vaf_log2r_plots/{sample_subfile}_somatic_varscan_filtered_NormADmax1-TumorADmin7-NormDPmin20_log2r.vcf", sample_subfile=config[sample]["tumors"]) for sample in SAMPLES]
        # [expand("vaf_log2r_plots/{sample_subfile}_vs_{normal}_somatic_varscan_filtered_NormADmax1-TumorADmin7-NormDPmin20_log2r.vcf", sample_subfile=config[sample]["tumors"], normal=config[sample]["normal"]) for sample in SAMPLES]
        [expand("vaf_log2r_plots/private_variants/plot_tables/{sample_subfile}_private_varscan_NormADmax1-NormDPmin20-TumorADmin7_log2r.txt", sample_subfile=config[sample]["tumors"]) for sample in SAMPLES]




#################################
##### Filter joint vcf file #####
#################################

#### Filter vcf file on number of variant reads
#### "FORMAT/AD[0:0]" refers to the AD value in the format column, first sample, first value,
#### which corresponds to the AD value of the variant allele in the blood sample

rule filter_joint_vcf:
    input:
        "{sample}_cns_varscan2.vcf.gz"
    output:
        "vaf_log2r_plots/{sample}_varscan_NormADmax1-NormDPmin20.vcf.gz"
    shell:
        """
        bcftools filter \
        --include 'FORMAT/AD[0:0] <= 1 && FILTER="PASS" && FORMAT/DP[0:0] >= 20' \
        {input} \
        -o {output}
        """


##################################################
##### Split joint vcf in one file per sample #####
##################################################

rule split_vcf:
    input:
        joint_vcf = lambda wildcards: expand("vaf_log2r_plots/{sample}_varscan_NormADmax1-NormDPmin20.vcf.gz", sample = '-'.join(wildcards.sample_subfile.split("-")[:2]))
    output:
        single_vcf = "vaf_log2r_plots/private_variants/{sample_subfile}_private_varscan_NormADmax1-NormDPmin20.vcf.gz"
    shell:
        """
        bcftools view  \
        -s {wildcards.sample_subfile}_tagseq-medexome \
        -o {output} \
        {input}
        """



###############################################
##### Prepare0 Log2R-data from ngcgh file #####
###############################################

#### Sort, bgzip and index ngcgh-file (required by bcftools)

rule sort_ngcgh:
    input:
        ngcgh_path + "{sample_subfile}_vs_{normal}-blod_normal_tagseq-medexome_ngcgh_w1000.txt"
    output:
        temp("vaf_log2r_plots/{sample_subfile}_vs_{normal}-blod_normal_tagseq-medexome_ngcgh_w1000_sorted.txt")
    shell:
        """
        sort -k1,1 -k2,2n -k3,3n {input} > {output}
        """

rule make_header:
    input:
        "vaf_log2r_plots/{sample_subfile}_vs_{normal}-blod_normal_tagseq-medexome_ngcgh_w1000_sorted.txt"
    output:
        temp("vaf_log2r_plots/{sample_subfile}_vs_{normal}-blod_normal_tagseq-medexome_ngcgh_w1000_sorted_header.txt")
    shell:
        """
        echo -e "CHROM\tFROM\tTO\tNORMCOUNT\tTUMORCOUNT\tLog2R" | \
        cat - {input} > {output}
        """

rule bgzip_ngcgh:
    input:
        "vaf_log2r_plots/{sample_subfile}_vs_{normal}-blod_normal_tagseq-medexome_ngcgh_w1000_sorted_header.txt"
    output:
        "vaf_log2r_plots/{sample_subfile}_vs_{normal}-blod_normal_tagseq-medexome_ngcgh_w1000_sorted_header.txt.gz"
    shell:
        """
        bgzip {input};
        """

rule index_ngcgh:
    input:
        "vaf_log2r_plots/{sample_subfile}_vs_{normal}-blod_normal_tagseq-medexome_ngcgh_w1000_sorted_header.txt.gz"
    output:
        "vaf_log2r_plots/{sample_subfile}_vs_{normal}-blod_normal_tagseq-medexome_ngcgh_w1000_sorted_header.txt.gz.tbi"
    shell:
        """
        tabix -s1 -b2 -e3 -S1 {input}
        """


###########################################
##### Annotate single-vcfs with log2R #####
###########################################

rule add_log2r_to_vcf:
    input:
        # ngcgh = expand("vaf_log2r_plots/{{sample_subfile}}_vs_{normal}-blod_normal_tagseq-medexome_ngcgh_w1000_sorted_header.txt.gz", normal=SAMPLES),
        # index = expand("vaf_log2r_plots/{{sample_subfile}}_vs_{normal}-blod_normal_tagseq-medexome_ngcgh_w1000_sorted_header.txt.gz.tbi", normal=SAMPLES),
        ngcgh = lambda wildcards: expand("vaf_log2r_plots/{{sample_subfile}}_vs_{normal}_tagseq-medexome_ngcgh_w1000_sorted_header.txt.gz", normal=config['-'.join(wildcards.sample_subfile.split("-")[:2])]["normal"]),
        index = lambda wildcards: expand("vaf_log2r_plots/{{sample_subfile}}_vs_{normal}_tagseq-medexome_ngcgh_w1000_sorted_header.txt.gz.tbi", normal=config['-'.join(wildcards.sample_subfile.split("-")[:2])]["normal"]),
        single_vcf = "vaf_log2r_plots/private_variants/{sample_subfile}_private_varscan_NormADmax1-NormDPmin20.vcf.gz"
    output:
        "vaf_log2r_plots/private_variants/{sample_subfile}_private_varscan_NormADmax1-NormDPmin20_log2r.vcf.gz"
    shell:
        """
        bcftools annotate \
        -a {input.ngcgh} \
        -c CHROM,FROM,TO,-,-,LOG2R \
        -h {header} \
        {input.single_vcf} \
        -o {output}
        """


###################################
##### Filter single vcf files #####
###################################

#### Filter vcf file on number of variant reads
#### "FORMAT/AD[0:0]" refers to the AD value in the format column, first sample, first value,
#### which corresponds to the AD value of the variant allele in the blood sample

rule filter_single_vcf:
    input:
        "vaf_log2r_plots/private_variants/{sample_subfile}_private_varscan_NormADmax1-NormDPmin20_log2r.vcf.gz"
    output:
        "vaf_log2r_plots/private_variants/{sample_subfile}_private_varscan_NormADmax1-NormDPmin20-TumorADmin7_log2r.vcf.gz"
    shell:
        """
        bcftools filter \
        --include 'FORMAT/AD[0:0] >= 7' \
        {input} \
        -o {output}
        """




############################################
##### Create table for VAF-Log2R plots #####
############################################

rule create_plot_table:
    input:
        "vaf_log2r_plots/private_variants/{sample_subfile}_private_varscan_NormADmax1-NormDPmin20-TumorADmin7_log2r.vcf.gz"
    output:
        "vaf_log2r_plots/private_variants/plot_tables/{sample_subfile}_private_varscan_NormADmax1-NormDPmin20-TumorADmin7_log2r.txt"
    shell:
        """
        bcftools query \
        -H \
        -f '%CHROM\t%POS\t[%FREQ\t]%INFO/LOG2R\t]\n' \
        {input} \
        -o {output}
        """
