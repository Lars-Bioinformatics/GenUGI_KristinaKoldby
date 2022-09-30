
SAMPLES, = glob_wildcards("{sample}-smoove.genotyped.vcf.gz")

rule all:
    input:
        expand("{sample}-smoove.filtered.vcf.gz", sample=SAMPLES)

rule remove_germline:
    input:
        "{sample}-smoove.genotyped.vcf.gz"
    output:
        "{sample}-smoove.filtered.vcf.gz"
    shell:
        """
        bcftools filter \
        --include 'FORMAT/AO[1:0] = 0' \
        {input} \
        -o {output}
        """
### From Lumpy vcf header:
##FORMAT=<ID=AO,Number=A,Type=Integer,Description="Alternate allele observations, with partial observations recorded fractionally">
