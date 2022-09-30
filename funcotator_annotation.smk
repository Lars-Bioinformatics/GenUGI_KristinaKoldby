### Before running script run command:
### Conda activate gatk4_update
### Otherwise bcftools reheader command will not work (requires bcftools v1.10)

### Makes vcf files a somewhat bigger, so run script on filtered vcf files
### Easiest way is to create a sub-folder and link relevant vcf-files in there
### (to avoid annotating filtered vcf-files in other sub-folders)

# SAMPLES, = glob_wildcards("{sample}.vcf.gz")
# SAMPLES, = glob_wildcards("{sample}.vcf")
# SAMPLES, = glob_wildcards("{sample}_newheader.vcf")
SAMPLES, = glob_wildcards("{sample}_pass.vcf")
# SAMPLES = "ECV2-4-biopsi-H1_tumor_tagseq-medexome-deep-seq"

ADmin = 3

fasta_fai = "/work/sdukoldby/resources/hg38/Homo_sapiens_assembly38.fasta.fai"
ref_build = "hg38"
ref_path = "/work/sdukoldby/resources/hg38/Homo_sapiens_assembly38.fasta"
data_source_path = "/work/sdukoldby/resources/exome_annotation"

### dbSNP-vcf-filen for hg38 har kromosomer angivet på formatet '1' i stedet for 'chr1' som i alle andre database-filer og variant-filter
### Har derfor reformateret dbSNPfilen "hg38_All_20170710.vcf.gz" til dette format og erstattet filen i data source mappen
### Command:
### bcftools view hg38_All_20170710.vcf.gz | awk '{if($0 !~ /^#/) print "chr"$0; else print $0}' > newfile
### herefter zipped med bgzip og indexeret med tabix

rule all:
    input:
        # expand("{sample}_funcotated_ADmin{ADmin}.txt", sample=SAMPLES, ADmin=ADmin)
        expand("{sample}_pass_funcotated_ADmin{ADmin}.txt", sample=SAMPLES, ADmin=ADmin)
        # expand("{sample}_pass_funcotated.txt", sample=SAMPLES)
        # expand("{sample}_funcotated.vcf", sample=SAMPLES)
        # expand("{sample}_newheader.vcf", sample=SAMPLES)

## Rename varscan vcf file - it's named .gz but it's not gzipped
# rule rename_varscan_vcf:
#     input:
#         gz = "{sample}.vcf.gz"
#     output:
#         vcf = "{sample}.vcf"
#     shell:
#         """
#         mv {input} {output}
#         """
#
# rule unzip_vcf:
#     input:
#         gz = "{sample}.vcf.gz"
#     output:
#         vcf = "{sample}.vcf"
#     shell:
#         """
#         bgzip -d {input}
#         """
#

### Prep varscan vcf
### Get chromosome length from reference file and add to header
# rule reformat_varscan_header:
#     input:
#         "{sample}.vcf"
#     output:
#         "{sample}_newheader.vcf"
#     shell:
#         """
#         bcftools reheader \
#         -f {fasta_fai} \
#         -o  {output}\
#         {input}
#         """
#
rule funcotator:
    input:
        # vcf = "{sample}_newheader.vcf"
        vcf = "{sample}_pass.vcf"
    output:
        vcf = "{sample}_pass_funcotated.vcf"
    shell:
        """
        gatk Funcotator \
        --variant {input} \
        --reference {ref_path} \
        --ref-version {ref_build} \
        --data-sources-path {data_source_path}  \
        --output {output} \
        --output-file-format VCF
        """

rule create_table:
    input:
        "{sample}_pass_funcotated.vcf"
    output:
        # "{sample}_funcotated_ADmin{ADmin}.txt"
        "{sample}_pass_funcotated_ADmin{ADmin}.txt"
    shell:
        """
        bcftools view \
        -i 'FORMAT/AD[*:1] >= {ADmin}' \
        {input} -Ou | \
        bcftools query \
        -H \
        -f '%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t%INFO/FUNCOTATION\t[%AF\t][%AD{{1}}\t][%DP\t]\n' \
        -o {output} \
        """

## AD{1} refers to the second value in the AD field, i.e. the variant allele

# For varscan files:
# shell:
#     """
#     bcftools view \
#     -i 'FORMAT/AD[*] >= {ADmin}' \
#     {input} -Ou | \
#     bcftools query \
#     -H \
#     -f '%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t%INFO/FUNCOTATION\t[%FREQ\t][%AD\t][%DP\t]\n' |\
#     sed -e 's/\%//g' > {output}
#     """
