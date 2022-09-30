# Rules in Snakefile based on https://github.com/cancerit/dockstore-cgpwgs/blob/develop/scripts/analysisWGS.sh

WORK="/work/sdukoldby/data/G45-2016_genugi/190325_A00653_0016_AHGL2LDSXX/BaseCalls/hg38"
INPUT_BAM = WORK+"/cgp_bam"
OUTPUT = WORK+"/wgs_out_hg38"

# SAMPLES, = glob_wildcards(INPUT+"{sample}_R1_001.fastq.gz")

# NORMAL = "G45-ECV2-31-blod_truseq-nano-genome_HGL2LDSXX_S27"
# TUMOR = "G45-ECV2-31-biopsi-F1_truseq-nano-genome_HGL2LDSXX_S22"

# NORMAL = "G56-blod-pt1_truseq-nano-genome_HGL2LDSXX_S5"
# TUMOR = "G56-sampleA1_truseq-nano-genome_HGL2LDSXX_S1"

configfile: WORK+"/samples_full_names.yaml"


# Resources - paths inside docker
ref = WORK+"/reference_files_GRCh38/core_ref_GRCh38_hla_decoy_ebv"
res = WORK+"/reference_files_GRCh38"

# Global variables
SPECIES = "Human"
ASSEMBLY = "NCBI38" #"NCBI37"
PROTOCOL = "WGS"
CAVESPLIT = 350000

# CONTIG_EXCLUDE = ""
CONTIG_EXCLUDE = "HLA%,chrUn%" # GRCh38: HLA%,chrUn%, GRCh37: NC_007605,hs37d5,GL%
if CONTIG_EXCLUDE == "":
    CAVE_CONTIG_EXCLUDE = ""
    PINDEL_CONTIG_EXCLUDE = ""
else:
    CAVE_CONTIG_EXCLUDE = "-x " + CONTIG_EXCLUDE
    PINDEL_CONTIG_EXCLUDE = "-e " + CONTIG_EXCLUDE


ASCAT_ADD_ARGS='' # Let ASCAT compute purity and ploidy
# ASCAT_ADD_ARGS='-pu ASCAT_PURITY -pi ASCAT_PLOIDY'
# ASCAT_ADD_ARGS='-pu 0.8 -pi 3'


onstart:
    shell("mkdir -p " + OUTPUT)
    shell("mkdir -p " + OUTPUT + "/finished")
    # shell("rm -f " + OUTPUT+"/finished/TEST_SUCCESS_"+TUMOR+"_vs_"+NORMAL+".txt")


rule all:
    input:
        [expand(OUTPUT+"/finished/SUCCESS_{tumor}_vs_{normal}.txt", tumor=config[s]["tumor"], normal=config[s]["normal"]) for s in config]
        # [expand(OUTPUT+"/finished/TEST_SUCCESS_{tumor}_vs_{normal}.txt", tumor=config[s]["tumor"], normal=config[s]["normal"]) for s in config]
        # expand(OUTPUT+"/finished/TEST_SUCCESS_{tumor}_vs_{normal}.txt", tumor=TUMOR, normal=NORMAL)
        # expand(OUTPUT+"/finished/SUCCESS_{tumor}_vs_{normal}.txt", tumor=TUMOR, normal=NORMAL)
        # expand(OUTPUT+"/{tumor}/timings/WGS_{tumor}_vs_{normal}.time.CaVEMan_setup", tumor=TUMOR, normal=NORMAL)


###############################################################################
#### Helper functions
###############################################################################
#def get_normalContamination(samplestats):
    #with open(samplestats) as f:
def get_normalContamination(wildcards):
    with open(checkpoints.ascat.get(**wildcards).output[1], 'r') as f:
        for line in f:
            if line.startswith("rho"):
                rho = line.strip().split(" ")
                # print(rho)
                break
    return(1-float(rho[1]))


###############################################################################
#### Parallel Block 1
#### CaVEMan setup
#### Genotype Check
#### VerifyBam Normal
###############################################################################
rule caveman_setup:
    input:
        normal=INPUT_BAM+"/{normal}/{normal}.bam",
        tumor=INPUT_BAM+"/{tumor}/{tumor}.bam"
    output:
        res1_1=OUTPUT+"/{tumor}/timings/WGS_{tumor}_vs_{normal}.time.CaVEMan_setup"
    threads: 24
    shell:
        """
        echo -e "[Parallel block 1] CaVEMan setup"; \
        TMP={OUTPUT}/{wildcards.tumor}/tmp; \
        mkdir -p {OUTPUT}/{wildcards.tumor}; \
        mkdir -p $TMP; \
    
        cp {res}/caveman/flag.vcf.config.WGS.ini {OUTPUT}/{wildcards.tumor}/flag.vcf.config.WGS.ini;
        perl -alne 'print join(qq{{\\t}},$F[0],0,$F[1],2);' < {ref}/genome.fa.fai | tee $TMP/norm.cn.bed > $TMP/tum.cn.bed;
        
        python2 ~/udocker-1.1.4/udocker run \
        --user=pi \
        --volume={WORK}:{WORK} \
        cgpwgs \
        bash -c '/usr/bin/time -v caveman.pl \
        -r {ref}/genome.fa.fai \
        -ig {res}/caveman/HiDepth.tsv \
        -b {res}/caveman/flagging \
        -ab {res}/vagrent \
        -u {res}/caveman \
        -s {SPECIES} \
        -sa {ASSEMBLY} \
        -t {threads} \
        -st {PROTOCOL} \
        -tc {OUTPUT}/{wildcards.tumor}/tmp/tum.cn.bed \
        -nc {OUTPUT}/{wildcards.tumor}/tmp/norm.cn.bed \
        -td 5 -nd 2 \
        -tb {input.tumor} \
        -nb {input.normal} \
        -c {OUTPUT}/{wildcards.tumor}/flag.vcf.config.WGS.ini \
        -f {res}/caveman/flagging/flag.to.vcf.convert.ini \
        -e {CAVESPLIT} \
        -o {OUTPUT}/{wildcards.tumor}/{PROTOCOL}_{wildcards.tumor}/caveman \
        {CAVE_CONTIG_EXCLUDE} \
        -p setup' \
        2>&1 | tee {output}
        """
        # 2>&1 | tee {output}
        # >& {output}
        # -x {CONTIG_EXCLUDE} \


rule genotype_check:
    input:
        normal=INPUT_BAM+"/{normal}/{normal}.bam",
        tumor=INPUT_BAM+"/{tumor}/{tumor}.bam"
    output:
        res1_2=OUTPUT+"/{tumor}/timings/WGS_{tumor}_vs_{normal}.time.geno"
    shell:
        """
        echo -e "[Parallel block 1] Genotype Check"; \
        python2 ~/udocker-1.1.4/udocker run \
        --user=laran \
        --volume={WORK}:{WORK} \
        cgpwgs \
        bash -c '/usr/bin/time -v compareBamGenotypes.pl \
        -o {OUTPUT}/{wildcards.tumor}/{PROTOCOL}_{wildcards.tumor}/genotyped \
        -j {OUTPUT}/{wildcards.tumor}/{PROTOCOL}_{wildcards.tumor}/genotyped/result.json \
        -nb {input.normal} \
        -tb {input.tumor} \
        -s {res}/qcGenotype/general.tsv \
        -g {res}/qcGenotype/gender.tsv' \
        2>&1 | tee {output}
        """


rule verifybam_normal:
    input:
        normal=INPUT_BAM+"/{normal}/{normal}.bam"
    output:
        res1_3=OUTPUT+"/{tumor}/timings/WGS_{normal}.time.verify_WT"
    threads: 24
    shell:
        """
        echo -e "[Parallel block 1] VerifyBam Normal"; \
        python2 ~/udocker-1.1.4/udocker run \
        --user=pi \
        --volume={WORK}:{WORK} \
        cgpwgs \
        bash -c '/usr/bin/time -v verifyBamHomChk.pl -d 25 \
        -o {OUTPUT}/{wildcards.tumor}/{PROTOCOL}_{wildcards.normal}/contamination \
        -b {input.normal} \
        -t {threads} \
        -j {OUTPUT}/{wildcards.tumor}/{PROTOCOL}_{wildcards.normal}/contamination/result.json \
        -s {res}/qcGenotype/verifyBamID_snps.vcf.gz' \
        2>&1 | tee {output}
        """
        # mv {OUTPUT}/{wildcards.tumor}/timings/{PROTOCOL}_{wildcards.tumor}.time.verify_WT {OUTPUT}/{wildcards.tumor}/timings/{PROTOCOL}_{wildcards.normal}.time.verify_WT


###############################################################################
#### Parallel Block 2
#### CaVEMan split
###############################################################################
rule caveman_split:
    input:
        normal=INPUT_BAM+"/{normal}/{normal}.bam",
        tumor=INPUT_BAM+"/{tumor}/{tumor}.bam",
        res1_1=OUTPUT+"/{tumor}/timings/WGS_{tumor}_vs_{normal}.time.CaVEMan_setup",
        res1_2=OUTPUT+"/{tumor}/timings/WGS_{tumor}_vs_{normal}.time.geno",
        res1_3=OUTPUT+"/{tumor}/timings/WGS_{normal}.time.verify_WT"
    output:
        res2_1=OUTPUT+"/{tumor}/timings/WGS_{tumor}_vs_{normal}.time.CaVEMan_split"
    threads: 24
    shell:
        """
        echo -e "[Parallel block 2] CaVEMan split"; \
        python2 ~/udocker-1.1.4/udocker run \
        --user=pi \
        --volume={WORK}:{WORK} \
        cgpwgs \
        bash -c '/usr/bin/time -v caveman.pl \
        -r {ref}/genome.fa.fai \
        -ig {res}/caveman/HiDepth.tsv \
        -b {res}/caveman/flagging \
        -ab {res}/vagrent \
        -u {res}/caveman \
        -s {SPECIES} \
        -sa {ASSEMBLY} \
        -t {threads} \
        -st {PROTOCOL} \
        -tc {OUTPUT}/{wildcards.tumor}/tmp/tum.cn.bed \
        -nc {OUTPUT}/{wildcards.tumor}/tmp/norm.cn.bed \
        -td 5 -nd 2 \
        -tb {input.tumor} \
        -nb {input.normal} \
        -c {OUTPUT}/{wildcards.tumor}/flag.vcf.config.WGS.ini \
        -f {res}/caveman/flagging/flag.to.vcf.convert.ini \
        -e {CAVESPLIT} \
        -o {OUTPUT}/{wildcards.tumor}/{PROTOCOL}_{wildcards.tumor}/caveman \
        {CAVE_CONTIG_EXCLUDE} \
        -p split' \
        2>&1 | tee {output}
        """


###############################################################################
#### Parallel Block 3
#### ASCAT
#### BRASS_input
#### BRASS_cover
###############################################################################
checkpoint ascat:
    input:
        normal=INPUT_BAM+"/{normal}/{normal}.bam",
        tumor=INPUT_BAM+"/{tumor}/{tumor}.bam",
        res2_1=OUTPUT+"/{tumor}/timings/WGS_{tumor}_vs_{normal}.time.CaVEMan_split"
    output:
        res3_1=OUTPUT+"/{tumor}/timings/WGS_{tumor}_vs_{normal}.time.ascat",
        samplestats=OUTPUT+"/{tumor}/WGS_{tumor}/ascat/{tumor}_vs_{normal}.samplestatistics.txt"
    # params:
    #     samplestats=OUTPUT+"/{tumor}/WGS_{tumor}/ascat/{tumor}.samplestatistics.txt"
    threads: 24
    shell:
        """
        echo -e "[Parallel block 3] ASCAT"; \
        python2 ~/udocker-1.1.4/udocker run \
        --user=pi \
        --volume={WORK}:{WORK} \
        cgpwgs \
        bash -c '/usr/bin/time -v ascat.pl \
        -o {OUTPUT}/{wildcards.tumor}/{PROTOCOL}_{wildcards.tumor}/ascat \
        -t {input.tumor} \
        -n {input.normal} \
        -sg {res}/ascat/SnpGcCorrections.tsv \
        -r {ref}/genome.fa \
        -q 20 \
        -g L \
        -l {res}/qcGenotype/gender.tsv \
        -rs {SPECIES} \
        -ra {ASSEMBLY} \
        -pr {PROTOCOL} \
        -pl ILLUMINA \
        -c {threads} \
        -force \
        {ASCAT_ADD_ARGS}' \
        2>&1 | tee {output.res3_1};

        cp {OUTPUT}/{wildcards.tumor}/WGS_{wildcards.tumor}/ascat/{wildcards.tumor}.samplestatistics.txt \
        {OUTPUT}/{wildcards.tumor}/WGS_{wildcards.tumor}/ascat/{wildcards.tumor}_vs_{wildcards.normal}.samplestatistics.txt
    
        cut -f 2,3,4,5 -d "," \
        {OUTPUT}/{wildcards.tumor}/WGS_{wildcards.tumor}/ascat/{wildcards.tumor}.copynumber.caveman.csv \
        | tr -s "," "\\t" > {OUTPUT}/{wildcards.tumor}/tmp/norm.cn.bed;
        cut -f 2,3,4,7 -d "," \
        {OUTPUT}/{wildcards.tumor}/WGS_{wildcards.tumor}/ascat/{wildcards.tumor}.copynumber.caveman.csv \
        | tr -s "," "\\t" > {OUTPUT}/{wildcards.tumor}/tmp/tum.cn.bed

        cp {OUTPUT}/{wildcards.tumor}/tmp/norm.cn.bed {OUTPUT}/{wildcards.tumor}/tmp/norm_ascat.cn.bed
        cp {OUTPUT}/{wildcards.tumor}/tmp/tum.cn.bed {OUTPUT}/{wildcards.tumor}/tmp/tum_ascat.cn.bed
        """
        # >& {output.res3_1};
        
        # -g XY \
        # -gc chrY \


rule BRASS_input:
    input:
        normal=INPUT_BAM+"/{normal}/{normal}.bam",
        tumor=INPUT_BAM+"/{tumor}/{tumor}.bam"
    output:
        res3_2=OUTPUT+"/{tumor}/timings/WGS_{tumor}_vs_{normal}.time.BRASS_input"
    threads: 24
    shell:
        """
        echo -e "[Parallel block 3] BRASS_input"; \
        python2 ~/udocker-1.1.4/udocker run \
        --user=pi \
        --volume={WORK}:{WORK} \
        cgpwgs \
        bash -c '/usr/bin/time -v brass.pl -j 4 -k 4 -c {threads} \
        -d {res}/brass/HiDepth.bed.gz \
        -f {res}/brass/brass_np.groups.gz \
        -g {ref}/genome.fa \
        -s {SPECIES} -as {ASSEMBLY} -pr {PROTOCOL} -pl ILLUMINA \
        -g_cache {res}/vagrent/vagrent.cache.gz \
        -vi {res}/brass/viral.genomic.fa.2bit \
        -mi {res}/brass/all_ncbi_bacteria \
        -b {res}/brass/500bp_windows.gc.bed.gz \
        -ct {res}/brass/CentTelo.tsv \
        -cb {res}/brass/cytoband.txt \
        -t {input.tumor} \
        -n {input.normal} \
        -o {OUTPUT}/{wildcards.tumor}/{PROTOCOL}_{wildcards.tumor}/brass \
        -p input'
        2>&1 | tee {output}
        """


rule BRASS_cover:
    input:
        normal=INPUT_BAM+"/{normal}/{normal}.bam",
        tumor=INPUT_BAM+"/{tumor}/{tumor}.bam"
    output:
        res3_3=OUTPUT+"/{tumor}/timings/WGS_{tumor}_vs_{normal}.time.BRASS_cover"
    threads: 24
    shell:
        """
        echo -e "[Parallel block 3] BRASS_cover"; \
        python2 ~/udocker-1.1.4/udocker run \
        --user=pi \
        --volume={WORK}:{WORK} \
        cgpwgs \
        bash -c '/usr/bin/time -v nice -n 10 brass.pl -j 4 -k 4 -c {threads} \
        -d {res}/brass/HiDepth.bed.gz \
        -f {res}/brass/brass_np.groups.gz \
        -g {ref}/genome.fa \
        -s {SPECIES} -as {ASSEMBLY} -pr {PROTOCOL} -pl ILLUMINA \
        -g_cache {res}/vagrent/vagrent.cache.gz \
        -vi {res}/brass/viral.genomic.fa.2bit \
        -mi {res}/brass/all_ncbi_bacteria \
        -b {res}/brass/500bp_windows.gc.bed.gz \
        -ct {res}/brass/CentTelo.tsv \
        -cb {res}/brass/cytoband.txt \
        -t {input.tumor} \
        -n {input.normal} \
        -o {OUTPUT}/{wildcards.tumor}/{PROTOCOL}_{wildcards.tumor}/brass \
        -p cover' \
        2>&1 | tee {output}
        """


###############################################################################
#### Parallel Block 4
#### cgpPindel
#### CaVEMan
###############################################################################
rule pindel:
    input:
        normal=INPUT_BAM+"/{normal}/{normal}.bam",
        tumor=INPUT_BAM+"/{tumor}/{tumor}.bam"
    output:
        res4_1=OUTPUT+"/{tumor}/timings/WGS_{tumor}_vs_{normal}.time.cgpPindel",
        bed=OUTPUT+"/{tumor}/WGS_{tumor}/pindel/{tumor}_vs_{normal}.germline.bed",
    threads: 24
    shell:
        """
        echo -e "[Parallel block 4] cgpPindel"; \
        python2 ~/udocker-1.1.4/udocker run \
        --user=pi \
        --volume={WORK}:{WORK} \
        cgpwgs \
        bash -c '/usr/bin/time -v pindel.pl \
        -o {OUTPUT}/{wildcards.tumor}/{PROTOCOL}_{wildcards.tumor}/pindel \
        -r {ref}/genome.fa \
        -t {input.tumor} \
        -n {input.normal} \
        -s {res}/pindel/simpleRepeats.bed.gz \
        -u {res}/pindel/pindel_np.gff3.gz \
        -f {res}/pindel/{PROTOCOL}_Rules.lst \
        -g {res}/vagrent/codingexon_regions.indel.bed.gz \
        -st {PROTOCOL} \
        -as {ASSEMBLY} \
        -sp {SPECIES} \
        {PINDEL_CONTIG_EXCLUDE} \
        -b {res}/pindel/HiDepth.bed.gz \
        -c {threads} \
        -sf {res}/pindel/softRules.lst' \
        2>&1 | tee {output.res4_1}
        """

rule caveman:
    input:
        normal=INPUT_BAM+"/{normal}/{normal}.bam",
        tumor=INPUT_BAM+"/{tumor}/{tumor}.bam",
        res3_1=OUTPUT+"/{tumor}/timings/WGS_{tumor}_vs_{normal}.time.ascat",
    output:
        res4_2=OUTPUT+"/{tumor}/timings/WGS_{tumor}_vs_{normal}.time.caveman"
    params:
        #samplestats=OUTPUT+"/{tumor}/WGS_{tumor}/ascat/{tumor}.samplestatistics.txt"
        # NORM_CONTAM=lambda wildcards: get_normalContamination(OUTPUT+"/"+
        #     wildcards.tumor+"/WGS_"+wildcards.tumor+
        #     "/ascat/"+wildcards.tumor+".samplestatistics.txt")
        # NORM_CONTAM=lambda wildcards: get_normalContamination(expand(OUTPUT+"/{tumor}/WGS_{tumor}/ascat/{tumor}.samplestatistics.txt", tumor=wildcards.tumor))
        NORM_CONTAM=get_normalContamination
    threads: 24
    shell:
        """
        echo -e "[Parallel block 4] CaVEMan"; \
        python2 ~/udocker-1.1.4/udocker run \
        --user=pi \
        --volume={WORK}:{WORK} \
        cgpwgs \
        bash -c '/usr/bin/time -v caveman.pl \
        -r {ref}/genome.fa.fai \
        -ig {res}/caveman/HiDepth.tsv \
        -b {res}/caveman/flagging \
        -ab {res}/vagrent \
        -u {res}/caveman \
        -s {SPECIES} \
        -sa {ASSEMBLY} \
        -t {threads} \
        -st {PROTOCOL} \
        -tc {OUTPUT}/{wildcards.tumor}/tmp/tum.cn.bed \
        -nc {OUTPUT}/{wildcards.tumor}/tmp/norm.cn.bed \
        -td 5 -nd 2 \
        -tb {input.tumor} \
        -nb {input.normal} \
        -c {OUTPUT}/{wildcards.tumor}/flag.vcf.config.WGS.ini \
        -f {res}/caveman/flagging/flag.to.vcf.convert.ini \
        -e {CAVESPLIT} \
        -o {OUTPUT}/{wildcards.tumor}/{PROTOCOL}_{wildcards.tumor}/caveman \
        {CAVE_CONTIG_EXCLUDE} \
        -k {params.NORM_CONTAM} \
        -no-flagging -noclean' \
        2>&1 | tee {output}
        """
        # NORM_CONTAM=`perl -ne 'if(m/^rho\s(.+)\n/){{print 1-$1;}}' {params.samplestats}`; \
        # -k ${{NORM_CONTAM}} \

###############################################################################
#### Parallel Block 5
#### BRASS
#### Pindel_annot
#### cgpFlagCaVEMan
###############################################################################
rule brass:
    input:
        normal=INPUT_BAM+"/{normal}/{normal}.bam",
        tumor=INPUT_BAM+"/{tumor}/{tumor}.bam",
        res3_1=OUTPUT+"/{tumor}/timings/WGS_{tumor}_vs_{normal}.time.ascat",
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
        -g {ref}/genome.fa \
        -s {SPECIES} -as {ASSEMBLY} -pr {PROTOCOL} -pl ILLUMINA \
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
        2>&1 | tee {output}
        """

rule Pindel_annot:
    input:
        res4_1=OUTPUT+"/{tumor}/timings/WGS_{tumor}_vs_{normal}.time.cgpPindel"
    output:
        res5_2=OUTPUT+"/{tumor}/timings/WGS_{tumor}_vs_{normal}.time.cgpPindel_annot"
    threads: 24
    shell:
        """
        echo -e "[Parallel block 5] Pindel_annot"

        # ensure no annotated pindel
        rm -f {OUTPUT}/{wildcards.tumor}/{PROTOCOL}_{wildcards.tumor}/pindel/{wildcards.tumor}_vs_{wildcards.normal}.annot.vcf.gz*

        python2 ~/udocker-1.1.4/udocker run \
        --user=pi \
        --volume={WORK}:{WORK} \
        cgpwgs \
        bash -c '/usr/bin/time -v AnnotateVcf.pl -t \
        -c {res}/vagrent/vagrent.cache.gz \
        -i {OUTPUT}/{wildcards.tumor}/{PROTOCOL}_{wildcards.tumor}/pindel/{wildcards.tumor}_vs_{wildcards.normal}.flagged.vcf.gz \
        -o {OUTPUT}/{wildcards.tumor}/{PROTOCOL}_{wildcards.tumor}/pindel/{wildcards.tumor}_vs_{wildcards.normal}.annot.vcf' \
        2>&1 | tee {output}
        """

rule cgpFlagCaVEMan:
    input:
        normal=INPUT_BAM+"/{normal}/{normal}.bam",
        tumor=INPUT_BAM+"/{tumor}/{tumor}.bam",
        bed=OUTPUT+"/{tumor}/WGS_{tumor}/pindel/{tumor}_vs_{normal}.germline.bed",
        res4_1=OUTPUT+"/{tumor}/timings/WGS_{tumor}_vs_{normal}.time.cgpPindel",
        res4_2=OUTPUT+"/{tumor}/timings/WGS_{tumor}_vs_{normal}.time.caveman",
    output:
        res5_3=OUTPUT+"/{tumor}/timings/WGS_{tumor}_vs_{normal}.time.cgpFlagCaVEMan"
    params:
        # samplestats=OUTPUT+"/{tumor}/WGS_{tumor}/ascat/{tumor}.samplestatistics.txt",
        # NORM_CONTAM=lambda wildcards: get_normalContamination(OUTPUT+"/"+
        #     wildcards.tumor+"/WGS_"+wildcards.tumor+
        #     "/ascat/"+wildcards.tumor+".samplestatistics.txt"),
        # GERMLINE_BED=lambda wildcards: expand(OUTPUT+"/{tumor}/WGS_{tumor}/pindel/{tumor}_vs_{normal}.germline.bed",
        #     tumor=wildcards.tumor, normal=wildcards.normal)
        NORM_CONTAM=get_normalContamination
    threads: 24
    shell:
        """
        echo -e "[Parallel block 5] CaVEMan flag"; \
        
        # need to sort and index pindel germline
        sort -k1,1 -k2,2n -k3,3n {input.bed} | bgzip -c > {input.bed}.gz
        tabix -p bed {input.bed}.gz
        
        python2 ~/udocker-1.1.4/udocker run \
        --user=pi \
        --volume={WORK}:{WORK} \
        cgpwgs \
        bash -c '/usr/bin/time -v caveman.pl \
        -r {ref}/genome.fa.fai \
        -ig {res}/caveman/HiDepth.tsv \
        -b {res}/caveman/flagging \
        -ab {res}/vagrent \
        -u {res}/caveman \
        -s {SPECIES} \
        -sa {ASSEMBLY} \
        -t {threads} \
        -st {PROTOCOL} \
        -tc {OUTPUT}/{wildcards.tumor}/tmp/tum.cn.bed \
        -nc {OUTPUT}/{wildcards.tumor}/tmp/norm.cn.bed \
        -td 5 -nd 2 \
        -tb {input.tumor} \
        -nb {input.normal} \
        -c {OUTPUT}/{wildcards.tumor}/flag.vcf.config.WGS.ini \
        -f {res}/caveman/flagging/flag.to.vcf.convert.ini \
        -e {CAVESPLIT} \
        -o {OUTPUT}/{wildcards.tumor}/{PROTOCOL}_{wildcards.tumor}/caveman \
        {CAVE_CONTIG_EXCLUDE} \
        -k {params.NORM_CONTAM} \
        -in {input.bed}.gz \
        -p flag' \
        2>&1 | tee {output}
        """
        # NORM_CONTAM=`perl -ne 'if(m/^rho\s(.+)\n/){{print 1-$1;}}' {params.samplestats}`;
        # rm -f $GERMLINE_BED
        
        # sort -k1,1 -k2,2n -k3,3n {params.GERMLINE_BED} | bgzip -c > {params.GERMLINE_BED}.gz
        # tabix -p bed {params.GERMLINE_BED}.gz
        
        # -in {params.GERMLINE_BED}.gz \

###############################################################################
#### Parallel Block 6
#### CaVEMan_annot
#### VerifyBam tumor
###############################################################################
rule caveman_annot:
    input:
        res5_3=OUTPUT+"/{tumor}/timings/WGS_{tumor}_vs_{normal}.time.cgpFlagCaVEMan"
    output:
        res6_1=OUTPUT+"/{tumor}/timings/WGS_{tumor}_vs_{normal}.time.CaVEMan_annot"
    threads: 24
    shell:
        """
        echo -e "[Parallel block 6] CaVEMan_annot";

        # annotate caveman
        rm -f {OUTPUT}/{wildcards.tumor}/{PROTOCOL}_{wildcards.tumor}/caveman/{wildcards.tumor}_vs_{wildcards.normal}.annot.muts.vcf.gz*
    
        python2 ~/udocker-1.1.4/udocker run \
        --user=pi \
        --volume={WORK}:{WORK} \
        cgpwgs \
        bash -c '/usr/bin/time -v AnnotateVcf.pl -t -c {res}/vagrent/vagrent.cache.gz \
        -i {OUTPUT}/{wildcards.tumor}/{PROTOCOL}_{wildcards.tumor}/caveman/{wildcards.tumor}_vs_{wildcards.normal}.flagged.muts.vcf.gz \
        -o {OUTPUT}/{wildcards.tumor}/{PROTOCOL}_{wildcards.tumor}/caveman/{wildcards.tumor}_vs_{wildcards.normal}.annot.muts.vcf' \
        2>&1 | tee {output}
        """

rule verifyBam_tumor:
    input:
        tumor=INPUT_BAM+"/{tumor}/{tumor}.bam",
        res3_1=OUTPUT+"/{tumor}/timings/WGS_{tumor}_vs_{normal}.time.ascat"
    output:
        res6_2=OUTPUT+"/{tumor}/timings/WGS_{tumor}_vs_{normal}.time.verify_MT"
    threads: 24
    shell:
        """
        echo -e "[Parallel block 6] VerifyBam Tumour"; \
        python2 ~/udocker-1.1.4/udocker run \
        --user=pi \
        --volume={WORK}:{WORK} \
        cgpwgs \
        bash -c '/usr/bin/time -v verifyBamHomChk.pl -d 25 \
        -o {OUTPUT}/{wildcards.tumor}/{PROTOCOL}_{wildcards.tumor}/contamination \
        -b {input.tumor} \
        -t {threads} \
        -a {OUTPUT}/{wildcards.tumor}/{PROTOCOL}_{wildcards.tumor}/ascat/{wildcards.tumor}.copynumber.caveman.csv \
        -j {OUTPUT}/{wildcards.tumor}/{PROTOCOL}_{wildcards.tumor}/contamination/result.json \
        -s {res}/qcGenotype/verifyBamID_snps.vcf.gz' \
        2>&1 | tee {output}
        """
        # mv {OUTPUT}/{wildcards.tumor}/timings/{PROTOCOL}_{wildcards.tumor}.time.verify_MT {OUTPUT}/{wildcards.tumor}/timings/{PROTOCOL}_{wildcards.tumor}.time.verify_MT


###############################################################################
#### Execution of battenberg copynumber (Manual)
###############################################################################
rule battenberg:
	input:
		normal=INPUT_BAM+"/{normal}/{normal}.bam",
		tumor=INPUT_BAM+"/{tumor}/{tumor}.bam",
		res2_1=OUTPUT+"/{tumor}/timings/WGS_{tumor}_vs_{normal}.time.CaVEMan_split"
	output:
		battenberg=OUTPUT+"/{tumor}/timings/WGS_{tumor}_vs_{normal}.time.battenberg"
	threads: 24
	shell:
		"""
		python2 ~/udocker-1.1.4/udocker run \
		--user=pi  \
		--volume={WORK}:{WORK}  \
		cgpwgs201 \
		bash -c '/usr/bin/time -v battenberg.pl \
		-o {OUTPUT}/{wildcards.tumor}/{PROTOCOL}_{wildcards.tumor}/battenberg \
		-r {ref}/genome.fa.fai \
		-u {res}/battenberg/1000genomesloci \
		-e {res}/battenberg/impute/impute_info.txt \
		-c {res}/battenberg/probloci.txt \
		-ig {res}/battenberg/ignore_contigs.txt \
		-gc {res}/battenberg/battenberg_wgs_gc_correction_1000g \
        -ge L \
        -gl {res}/qcGenotype/gender.tsv \
		-tb {input.tumor} \
		-nb {input.normal} \
		-rs {SPECIES} \
		-ra {ASSEMBLY} \
		-nl 50 \
		-t {threads}' \
		2>&1 | tee {output}

		zcat {OUTPUT}/{wildcards.tumor}/{PROTOCOL}_{wildcards.tumor}/battenberg/{wildcards.tumor}_battenberg_cn.vcf.gz | \
		grep -v ^# | cut -f 1,2,8,10 | \
		awk "{{for(i=1;i<=NF;i++){{if (i==3) {{split($i,a,"=");printf("%s\\t",a[3])}} else if (i==4) {{split($i,a,":");printf("%s\\n", a[2])}} else {{printf("%s\\t",$i)}} }}}}" > \
		{OUTPUT}/{wildcards.tumor}/tmp/norm_bat.cn.bed

        zcat {OUTPUT}/{wildcards.tumor}/{PROTOCOL}_{wildcards.tumor}/battenberg/{wildcards.tumor}_battenberg_cn.vcf.gz | \
        grep -v ^# | cut -f 1,2,8,11 | \
        awk "{{for(i=1;i<=NF;i++){{if (i==3) {{split($i,a,"=");printf("%s\\t",a[3])}} else if (i==4) {{split($i,a,":");printf("%s\\n", a[2])}} else {{printf("%s\\t",$i)}} }}}}" > \
        {OUTPUT}/{wildcards.tumor}/tmp/tum_bat.cn.bed
		"""
        # -ge L \
        # -gl {res}/qcGenotype/gender.tsv \


###############################################################################
#### Tying it all together with max possible paralization
###############################################################################
rule run_parallel:
    input:
        res1_3=OUTPUT+"/{tumor}/timings/WGS_{normal}.time.verify_WT",
        res5_1=OUTPUT+"/{tumor}/timings/WGS_{tumor}_vs_{normal}.time.BRASS",
        res5_2=OUTPUT+"/{tumor}/timings/WGS_{tumor}_vs_{normal}.time.cgpPindel_annot",
        res6_1=OUTPUT+"/{tumor}/timings/WGS_{tumor}_vs_{normal}.time.CaVEMan_annot",
        res6_2=OUTPUT+"/{tumor}/timings/WGS_{tumor}_vs_{normal}.time.verify_MT",
        # battenberg=OUTPUT+"/{tumor}/timings/WGS_{tumor}_vs_{normal}.time.battenberg"
    output:
        OUTPUT+"/finished/SUCCESS_{tumor}_vs_{normal}.txt"
    shell:
        """
        echo 'Successfully completed entire cgp pipeline!' > {output}
        """

rule run_test:
    input:
        res1_1=OUTPUT+"/{tumor}/timings/WGS_{tumor}_vs_{normal}.time.CaVEMan_setup",
        # res1_2=OUTPUT+"/{tumor}/timings/WGS_{tumor}_vs_{normal}.time.geno",
        # res1_3=OUTPUT+"/{tumor}/timings/WGS_{normal}.time.verify_WT",
        # res5_1=OUTPUT+"/{tumor}/timings/WGS_{tumor}_vs_{normal}.time.BRASS",
        # res5_2=OUTPUT+"/{tumor}/timings/WGS_{tumor}_vs_{normal}.time.cgpPindel_annot",
        # res6_1=OUTPUT+"/{tumor}/timings/WGS_{tumor}_vs_{normal}.time.CaVEMan_annot",
        # res6_2=OUTPUT+"/{tumor}/timings/WGS_{tumor}_vs_{normal}.time.verify_MT",
    output:
        OUTPUT+"/finished/TEST_SUCCESS_{tumor}_vs_{normal}.txt"
    shell:
        """
        echo 'Successfully completed entire cgp pipeline!' > {output}
        """
    