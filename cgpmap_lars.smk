WORK="/work/sdukoldby/data/G45-2016_genugi/190325_A00653_0016_AHGL2LDSXX/BaseCalls/"
INPUT = WORK+"fastq_files_for_sanger/"
OUTPUT = WORK+"hg38/cgp_bam/"
YAML = WORK+"hg38/ReadGroupsYAML/"

SAMPLES, = glob_wildcards(INPUT+"{sample}_1.fastq.gz")
# SAMPLES = "G45-ECV2-29-biopsi-C1_truseq-nano-genome_HGL2LDSXX_S21"
# SAMPLES = "G45-ECV2-35-blod_truseq-nano-genome_HGL2LDSXX_S28"

# Resources - paths inside docker
ref = WORK+"hg38/reference_files_GRCh38/core_ref_GRCh38_hla_decoy_ebv"
# ref = WORK+"/reference_files_GRCh37/"

rule all:
    input:
        #expand(OUTPUT+"{sample}/{sample}.bam", sample=SAMPLES)
        expand(OUTPUT+"{sample}/{sample}.bam", sample=SAMPLES)

## IMPORTANT NOTE: Make sure fastq files ends on _1.fastq.gz and _2.fastq.gz
## - otherwise faulty bam files are produced
rule cgpmap:
    input:
        f1=INPUT+"{sample}_1.fastq.gz",
        f2=INPUT+"{sample}_2.fastq.gz",
        rg=YAML+"ReadGroups_{sample}.yaml"
    output:
        bam=OUTPUT+"{sample}/{sample}.bam"
    threads: 24
    shell:
        """
        python2 ~/udocker-1.1.4/udocker run \
        --user=laran \
        --volume={WORK}:{WORK} \
        cgpmap \
        ds-cgpmap.pl \
        -reference {ref} \
        -bwa_idx {ref} \
        -sample {wildcards.sample} \
        -groupinfo {input.rg} \
        -threads {threads} \
        -outdir {OUTPUT}{wildcards.sample} \
        {input.f1} {input.f2}
        """