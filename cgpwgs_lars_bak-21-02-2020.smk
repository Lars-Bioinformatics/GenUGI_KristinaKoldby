# Rules in Snakefile based on https://github.com/cancerit/dockstore-cgpwgs/blob/develop/scripts/analysisWGS.sh

WORK="/work/sdukoldby/data/G45-2016_genugi/190325_A00653_0016_AHGL2LDSXX/BaseCalls"
INPUT_BAM = WORK+"/cgp_bam"
OUTPUT = WORK+"/hg38/wgs_out_hg38"

# SAMPLES, = glob_wildcards(INPUT+"{sample}_R1_001.fastq.gz")
# SAMPLES = "G45-ECV2-29-biopsi-C1_truseq-nano-genome_HGL2LDSXX_S21"

# NORMAL = "G45-ECV2-31-blod_truseq-nano-genome_HGL2LDSXX_S27"
# TUMOR = "G45-ECV2-31-biopsi-F1_truseq-nano-genome_HGL2LDSXX_S22"
# NORMAL = "G45-ECV2-31-blod"
# TUMOR = "G45-ECV2-31-biopsi-F1"

# Resources - paths inside docker
# ref = WORK+"/reference_files_GRCh38/core_ref_GRCh38_hla_decoy_ebv"
# res = WORK+"/reference_files_GRCh38"
# ref = WORK+"/reference_files_hg38/core_ref_GRCh38_hla_decoy_ebv"
# ref = WORK+"/resources/hg38"
ref = WORK+"/hg38/reference_files_GRCh38/core_ref_GRCh38_hla_decoy_ebv"
res = WORK+"/hg38/reference_files_GRCh38"


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
        perl -alne 'print join(qq{{\\t}},$F[0],0,$F[1],2);' < {ref}/Homo_sapiens_assembly38.fasta.fai | tee $TMP/norm.cn.bed > $TMP/tum.cn.bed;
        
        python2 ~/udocker-1.1.4/udocker run \
        --user=pi \
        --volume={WORK}:{WORK} \
        cgpwgs200 \
        bash -c '/usr/bin/time -v caveman.pl \
        -r {ref}/Homo_sapiens_assembly38.fasta.fai \
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
        -x {CONTIG_EXCLUDE} \
        -p setup' \
        2>&1 | tee {output}
        """
        # >& {output}


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
        -s {res}/qcGenotype_GRCh38_hla_decoy_ebv/general.tsv \
        -g {res}/qcGenotype_GRCh38_hla_decoy_ebv/gender.tsv' \
        >& {output}
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
        -s {res}/qcGenotype_GRCh38_hla_decoy_ebv/verifyBamID_snps.vcf.gz' \
        >& {output}
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
        -r {ref}/Homo_sapiens_assembly38.fasta.fai \
        -ig {res}/caveman/HiDepth.tsv \
        -b {res}/caveman/flagging \
        -ab {res}/vagrent \
        -u {res}/caveman \
        -s \'{SPECIES}\' \
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
        -x {CONTIG_EXCLUDE} \
        -p split' \
        >& {output}
        """


###############################################################################
#### Parallel Block 3
#### ASCAT
#### BRASS_input
#### BRASS_cover
###############################################################################
rule ascat:
    input:
        normal=INPUT_BAM+"/{normal}/{normal}.bam",
        tumor=INPUT_BAM+"/{tumor}/{tumor}.bam",
        res2_1=OUTPUT+"/{tumor}/timings/WGS_{tumor}_vs_{normal}.time.CaVEMan_split"
    output:
        res3_1=OUTPUT+"/{tumor}/timings/WGS_{tumor}_vs_{normal}.time.ascat",
        samplestats=OUTPUT+"/{tumor}/WGS_{tumor}/ascat/{normal}.samplestatistics.txt"
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
        -r {ref}/Homo_sapiens_assembly38.fasta.fai \
        -q 20 \
        -g L \
        -l {res}/qcGenotype_GRCh38_hla_decoy_ebv/gender.tsv \
        -rs \'{SPECIES}\' \
        -ra {ASSEMBLY} \
        -pr {PROTOCOL} \
        -pl ILLUMINA \
        -c {threads} \
        -force \
        {ASCAT_ADD_ARGS}' \
        >& {output.res3_1};
    
        cut -f 2,3,4,5 -d "," \
        {OUTPUT}/{wildcards.tumor}/WGS_{wildcards.tumor}/ascat/{wildcards.tumor}.copynumber.caveman.csv \
        | tr -s "," "\\t" > /work/wgs_out/{wildcards.tumor}/tmp/norm_ascat.cn.bed;
        cut -f 2,3,4,7 -d "," \
        {OUTPUT}/{wildcards.tumor}/WGS_{wildcards.tumor}/ascat/{wildcards.tumor}.copynumber.caveman.csv \
        | tr -s "," "\\t" > /work/wgs_out/{wildcards.tumor}/tmp/tum_ascat.cn.bed
        """
        # >& {output.res3_1};
        
        # -g XY \
        # -gc chrY \
        
        # -g L \
        # -l {res}/qcGenotype_GRCh38_hla_decoy_ebv/gender.tsv \
        
        


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
        -o {OUTPUT}/{wildcards.tumor}/{PROTOCOL}_{wildcards.tumor}/brass \
        -p input'
        >& {output}
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
        -o {OUTPUT}/{wildcards.tumor}/{PROTOCOL}_{wildcards.tumor}/brass \
        -p cover' \
        >& {output}
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
        res4_1=OUTPUT+"/{tumor}/timings/WGS_{tumor}_vs_{normal}.time.cgpPindel"
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
        -r {ref}/Homo_sapiens_assembly38.fasta \
        -t {input.tumor} \
        -n {input.normal} \
        -s {res}/pindel/simpleRepeats.bed.gz \
        -u {res}/pindel/pindel_np.gff3.gz \
        -f {res}/pindel/{PROTOCOL}_Rules.lst \
        -g {res}/vagrent/codingexon_regions.indel.bed.gz \
        -st {PROTOCOL} \
        -as {ASSEMBLY} \
        -sp '{SPECIES}' \
        -e {CONTIG_EXCLUDE} \
        -b {res}/pindel/HiDepth.bed.gz \
        -c {threads} \
        -sf {res}/pindel/softRules.lst' \
        >& {output}
        """


rule caveman:
    input:
        normal=INPUT_BAM+"/{normal}/{normal}.bam",
        tumor=INPUT_BAM+"/{tumor}/{tumor}.bam",
        samplestats=OUTPUT+"/{tumor}/WGS_{tumor}/ascat/{normal}.samplestatistics.txt",
        #res3_1=OUTPUT+"/{tumor}/timings/WGS_{tumor}_vs_{normal}.time.ascat"
    output:
        res4_2=OUTPUT+"/{tumor}/timings/WGS_{tumor}_vs_{normal}.time.caveman"
    threads: 24
    shell:
        """
        echo -e "[Parallel block 4] CaVEMan"; \
        NORM_CONTAM=`perl -ne 'if(m/^rho\s(.+)\n/){{print 1-$1;}}' {input.samplestats}`; \
        python2 ~/udocker-1.1.4/udocker run \
        --user=pi \
        --volume={WORK}:{WORK} \
        cgpwgs \
        bash -c '/usr/bin/time -v caveman.pl \
        -r {ref}/Homo_sapiens_assembly38.fasta.fai \
        -ig {res}/caveman/HiDepth.tsv \
        -b {res}/caveman/flagging \
        -ab {res}/vagrent \
        -u {res}/caveman \
        -s '{SPECIES}' \
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
        -x {CONTIG_EXCLUDE} \
        -k $NORM_CONTAM \
        -no-flagging -noclean' \
        >& {output}
        """

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
        rm -f {OUTPUT}/{wildcards.tumor}/{PROTOCOL}_{wildcards.tumor}/pindel/{wildcards.tumor}.annot.vcf.gz*

        python2 ~/udocker-1.1.4/udocker run \
        --user=pi \
        --volume={WORK}:{WORK} \
        cgpwgs \
        bash -c '/usr/bin/time -v AnnotateVcf.pl -t \
        -c {res}/vagrent/vagrent.cache.gz \
        -i {OUTPUT}/{wildcards.tumor}/{PROTOCOL}_{wildcards.tumor}/pindel/{wildcards.tumor}.flagged.vcf.gz \
        -o {OUTPUT}/{wildcards.tumor}/{PROTOCOL}_{wildcards.tumor}/pindel/{wildcards.tumor}.annot.vcf' \
        >& {output}
        """

rule cgpFlagCaVEMan:
    input:
        normal=INPUT_BAM+"/{normal}/{normal}.bam",
        tumor=INPUT_BAM+"/{tumor}/{tumor}.bam",
        samplestats=OUTPUT+"/{tumor}/WGS_{tumor}/{tumor}//ascat/{normal}.samplestatistics.txt",
        res4_2=OUTPUT+"/{tumor}/timings/WGS_{tumor}_vs_{normal}.time.caveman"
    output:
        res5_3=OUTPUT+"/{tumor}/timings/WGS_{tumor}_vs_{normal}.time.cgpFlagCaVEMan"
    threads: 24
    shell:
        """
        echo -e "[Parallel block 5] CaVEMan flag"; \
    
        NORM_CONTAM=`perl -ne 'if(m/^rho\s(.+)\n/){{print 1-$1;}}' {input.samplestats}`;
        
        GERMLINE_BED={OUTPUT}/{wildcards.tumor}/{PROTOCOL}_{wildcards.tumor}/pindel/{wildcards.tumor}.germline.bed
        # need to sort and index pindel germline
        sort -k1,1 -k2,2n -k3,3n $GERMLINE_BED | bgzip -c > $GERMLINE_BED.gz
        tabix -p bed $GERMLINE_BED.gz
        rm -f $GERMLINE_BED
        
        python2 ~/udocker-1.1.4/udocker run \
        --user=pi \
        --volume={WORK}:{WORK} \
        cgpwgs \
        bash -c '/usr/bin/time -v caveman.pl \
        -r {ref}/Homo_sapiens_assembly38.fasta.fai \
        -ig {res}/caveman/HiDepth.tsv \
        -b {res}/caveman/flagging \
        -ab {res}/vagrent \
        -u {res}/caveman \
        -s '{SPECIES}' \
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
        -x {CONTIG_EXCLUDE} \
        -k $NORM_CONTAM \
        -in $GERMLINE_BED.gz \
        -p flag' \
        >& {output}
        """

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
        rm -f {OUTPUT}/{PROTOCOL}_{wildcards.tumor}/caveman/{wildcards.tumor}.annot.muts.vcf.gz*
    
        python2 ~/udocker-1.1.4/udocker run \
        --user=pi \
        --volume={WORK}:{WORK} \
        cgpwgs \
        bash -c '/usr/bin/time -v AnnotateVcf.pl -t -c {res}/vagrent/vagrent.cache.gz \
        -i {OUTPUT}/{PROTOCOL}_{wildcards.tumor}/caveman/{wildcards.tumor}.flagged.muts.vcf.gz \
        -o {OUTPUT}/{PROTOCOL}_{wildcards.tumor}/caveman/{wildcards.tumor}.annot.muts.vcf' \
        >& {output}
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
        -s {res}/verifyBamID_snps.vcf.gz' \
        >& {output}
        """
        # mv {OUTPUT}/{wildcards.tumor}/timings/{PROTOCOL}_{wildcards.tumor}.time.verify_MT {OUTPUT}/{wildcards.tumor}/timings/{PROTOCOL}_{wildcards.tumor}.time.verify_MT

###############################################################################
#### Tying it all together with max possible paralization
###############################################################################
rule run_parallel:
    input:
        res1_3=OUTPUT+"/{tumor}/timings/WGS_{normal}.time.verify_WT",
        res5_1=OUTPUT+"/{tumor}/timings/WGS_{tumor}_vs_{normal}.time.BRASS",
        res5_2=OUTPUT+"/{tumor}/timings/WGS_{tumor}_vs_{normal}.time.cgpPindel_annot",
        res6_1=OUTPUT+"/{tumor}/timings/WGS_{tumor}_vs_{normal}.time.CaVEMan_annot",
        res6_2=OUTPUT+"/{tumor}/timings/WGS_{tumor}.time.verify_MT",
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
    