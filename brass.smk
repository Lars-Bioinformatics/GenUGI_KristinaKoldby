# Rules in Snakefile based on https://github.com/cancerit/dockstore-cgpwgs/blob/develop/scripts/analysisWGS.sh

WORK="/work/sdukoldby"
INPUT_BAM = WORK+"/data/G45-2016_genugi/190325_A00653_0016_AHGL2LDSXX/BaseCalls/cgp_bam"
OUTPUT = WORK+"/data/G45-2016_genugi/190325_A00653_0016_AHGL2LDSXX/BaseCalls/hg38/wgs_out_hg38"

# SAMPLES, = glob_wildcards(INPUT+"{sample}_R1_001.fastq.gz")
# SAMPLES = "G45-ECV2-29-biopsi-C1_truseq-nano-genome_HGL2LDSXX_S21"

NORMAL = "G45-ECV2-31-blod_truseq-nano-genome_HGL2LDSXX_S27"
TUMOR = "G45-ECV2-31-biopsi-F1_truseq-nano-genome_HGL2LDSXX_S22"
# NORMAL = "G45-ECV2-31-blod"
# TUMOR = "G45-ECV2-31-biopsi-F1"


# Resources - paths inside docker
# ref = WORK+"/reference_files_GRCh38/core_ref_GRCh38_hla_decoy_ebv"
# res = WORK+"/reference_files_GRCh38"
# ref = WORK+"/reference_files_hg38/core_ref_GRCh38_hla_decoy_ebv"
ref = WORK+"/resources/hg38"
res = WORK+"/data/G45-2016_genugi/190325_A00653_0016_AHGL2LDSXX/BaseCalls/hg38/reference_files_GRCh38"


# Global variables
SPECIES = "Human"
ASSEMBLY = "NCBI38"
PROTOCOL = "WGS"
CAVESPLIT = 350000
CONTIG_EXCLUDE = "" # HLA*,chrUn*

ASCAT_ADD_ARGS='' # Let ASCAT compute purity and ploidy
# ASCAT_ADD_ARGS='-pu ASCAT_PURITY -pi ASCAT_PLOIDY'
# ASCAT_ADD_ARGS='-pu 0.8 -pi 3'


onstart:
    shell("mkdir -p " + OUTPUT)
    shell("mkdir -p " + OUTPUT + "/finished")
    shell("rm -f " + OUTPUT+"/finished/TEST_SUCCESS_"+TUMOR+"_vs_"+NORMAL+".txt")


rule all:
    input:
        expand(OUTPUT+"/finished/TEST_SUCCESS_{tumor}_vs_{normal}.txt", tumor=TUMOR, normal=NORMAL)
        # expand(OUTPUT+"/finished/SUCCESS_{tumor}_vs_{normal}.txt", tumor=TUMOR, normal=NORMAL)
        # expand(OUTPUT+"/{tumor}/timings/WGS_{tumor}_vs_{normal}.time.CaVEMan_setup", tumor=TUMOR, normal=NORMAL)

rule brass:
    input:
        normal=INPUT_BAM+"/{normal}/{normal}.bam",
        tumor=INPUT_BAM+"/{tumor}/{tumor}.bam",
        res3_2=OUTPUT+"/{tumor}/timings/WGS_{tumor}_vs_{normal}.time.BRASS_input",
        res3_3=OUTPUT+"/{tumor}/timings/WGS_{tumor}_vs_{normal}.time.BRASS_cover"
    output:
        res5_1=OUTPUT+"/{tumor}/timings/WGS_{tumor}_vs_{normal}.time.BRASS"
    threads: 24
    shell:
        """
        echo -e "[Parallel block 5] BRASS"; \
        python2 ~/udocker-1.1.4/udocker run \
        --user=pi \
        --volume={WORK}:{WORK} \
        cgpwgs \
        bash -c '/usr/bin/time -v brass.pl -j 4 -k 4 -c {threads} \
        -d {res}/brass/HiDepth.bed.gz \
        -f {res}/brass/brass_np.groups.gz \
        -g {ref}/Homo_sapiens_assembly38.fasta \
        -s '{SPECIES}' -as {ASSEMBLY} -pr {PROTOCOL} -pl ILLUMINA \
        -g_cache {res}/vagrent/vagrent.cache.gz \
        -vi {res}/brass/viral.genomic.fa.2bit \
        -mi {res}/brass/all_ncbi_bacteria \
        -b {res}/brass/500bp_windows.gc.bed.gz \
        -ct {res}/brass/CentTelo.tsv \
        -cb {res}/brass/cytoband.txt \
        -t {input.tumor} \
        -n {input.normal} \
        -ss {OUTPUT}/{wildcards.tumor}/{PROTOCOL}_{wildcards.tumor}/ascat/{wildcards.tumor}.samplestatistics.txt \
        -o {OUTPUT}/{wildcards.tumor}/{PROTOCOL}_{wildcards.tumor}/brass' \
        >& {output}
        """
