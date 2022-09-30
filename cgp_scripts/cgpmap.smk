INPUT = "/work/sdukoldby/data/G45-2016_genugi/190325_A00653_0016_AHGL2LDSXX/BaseCalls/"
WORK="/work/sdukoldby/data/G45-2016_genugi/wgs_bam/"

# SAMPLES, = glob_wildcards(INPUT+"{sample}_R1_001.fastq.gz")
SAMPLES = "G45-ECV2-29-biopsi-C1_truseq-nano-genome_HGL2LDSXX_S21"

# Resources - paths inside docker
# ref = "/work/reference_files_GRCh38/archives/core_ref_GRCh38_hla_decoy_ebv.tar.gz"
# bwa_index = "/work/reference_files_GRCh38/archives/bwa_idx_GRCh38_hla_decoy_ebv.tar.gz"
ref = WORK+"/reference_files_GRCh38/core_ref_GRCh38_hla_decoy_ebv"

# onstart:
#     shell("conda activate py2")

rule all:
    input:
        expand(WORK+"cgp_out/{sample}.bam", sample=SAMPLES)

rule cgpmap:
    input:
        INPUT+"{sample}_R1_001.fastq.gz",
        INPUT+"{sample}_R2_001.fastq.gz"
    output:
        WORK+"cgp_out/{sample}.bam"
    shell:
        """
        python2 ~/udocker-1.1.4/udocker run \
        --user=kmkoldby  \
        --volume={WORK}:{WORK} \
        --volume={INPUT}:{INPUT} \
        cgpmap \
        ds-cgpmap.pl \
        -r {ref} \
        -i {ref} \
        -s {wildcards.sample}  \
        -t 24 \
        -o {WORK}/cgp_out \
        {input}
        """
        #/work/fastq/{wildcards.sample}_R[12]_001.fastq.gz
        # --volume={WORK}:/work \
