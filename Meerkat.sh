#!/bin/sh
#
#SBATCH --account sdukoldby_slim      # account
#SBATCH --nodes 1                 # number of nodes
#SBATCH --time 24:00:00            # max time (HH:MM:SS)

### Echo on
set -x

#SAMPLEID=$1
SAMPLEID="G45-ECV2-4-blod_truseq-nano-genome_HGL2LDSXX_S26"
#SAMPLEID="G45-ECV2-8-blod_truseq-nano-genome_HGL2LDSXX_S29"
#SAMPLEID="G45-ECV2-29-blod_truseq-nano-genome_HGL2LDSXX_S30"
#SAMPLEID="G45-ECV2-31-blod_truseq-nano-genome_HGL2LDSXX_S27"
#SAMPLEID="G45-ECV2-35-blod_truseq-nano-genome_HGL2LDSXX_S28"
#SAMPLEID="G45-ECV2-4-biopsi-H2_truseq-nano-genome_HGL2LDSXX_S24"
#SAMPLEID="G45-ECV2-8-biopsi-I2_truseq-nano-genome_HGL2LDSXX_S25"
#SAMPLEID="G45-ECV2-29-biopsi-C1_truseq-nano-genome_HGL2LDSXX_S21"
#SAMPLEID="G45-ECV2-31-biopsi-F1_truseq-nano-genome_HGL2LDSXX_S22"
#SAMPLEID="G45-ECV2-35-biopsi-G1_truseq-nano-genome_HGL2LDSXX_S23"

TUMORID="G45-ECV2-4-biopsi-H2_truseq-nano-genome_HGL2LDSXX_S24"
#TUMORID="G45-ECV2-8-biopsi-I2_truseq-nano-genome_HGL2LDSXX_S25"
#TUMORID="G45-ECV2-29-biopsi-C1_truseq-nano-genome_HGL2LDSXX_S21"
#TUMORID="G45-ECV2-31-biopsi-F1_truseq-nano-genome_HGL2LDSXX_S22"
#TUMORID="G45-ECV2-35-biopsi-G1_truseq-nano-genome_HGL2LDSXX_S23"
TUMORID2="G45-ECV2-4-biopsi-H2"
#TUMORID2="G45-ECV2-8-biopsi-I2"
#TUMORID2="G45-ECV2-29-biopsi-C1"
#TUMORID2="G45-ECV2-31-biopsi-F1"
#TUMORID2="G45-ECV2-35-biopsi-G1"

MER="/work/sdukoldby/tools/Meerkat"
SCRIPTS="/work/sdukoldby/tools/Meerkat/scripts"
DAT="/work/sdukoldby/data/G45-2016_genugi/190325_A00653_0016_AHGL2LDSXX/BaseCalls/hg38/meerkat"
OUT="/work/sdukoldby/data/G45-2016_genugi/190325_A00653_0016_AHGL2LDSXX/BaseCalls/hg38/meerkat"
BLAST="/work/sdukoldby/tools/ncbi-blast-2.9.0+/bin"
REF="/work/sdukoldby/resources/hg38"
BIODB="/work/sdukoldby/tools/miniconda3/lib/site_perl/5.26.2/Bio/DB"

### Run Meerkat
### When running Meerkat on tumor samples: Run pre_process.pl as usual. Then replace tumor.blacklist.gz in data folder by blacklist.gz file generated for matched normal sample and rename this file using prefix for tumor sample.
###REMEMBER TO ACTIVATE CONDA ENVIRONMENT "meerkat" ###

export LD_LIBRARY_PATH=/work/sdukoldby/tools/Meerkat/src/mybamtools/lib
cd $DAT
#echo $SAMPLEID
#echo $PWD

#perl $SCRIPTS/pre_process.pl -b $SAMPLEID.recalibrated.bam -I $REF/Homo_sapiens_assembly38.fasta -A $REF/Homo_sapiens_assembly38.fasta.fai -t 24 -s 20 -k 1500 -q 15 -l 0
#perl $SCRIPTS/meerkat.pl -b $SAMPLEID.recalibrated.bam -F $REF -B $BLAST -t 24 -s 20 -m 0 -d 5 -c 5 -o 1 -p 3 -l 0
#perl $SCRIPTS/mechanism.pl -b $SAMPLEID.recalibrated.bam -R $REF/rmsk-hg38.txt

###Somatic filtering###
##### $SAMPLEID should be the normal sample matching the tumor sample ### 
perl $SCRIPTS/somatic_sv.pl -i $TUMORID/$TUMORID.recalibrated.variants -o $TUMORID/$TUMORID2.meerkat.somatic1.variants -F $DAT/normal_discord_files -l 1000 -R $REF/rmsk-hg38.txt
#perl $SCRIPTS/somatic_sv.pl -i $TUMORID/$TUMORID2.meerkat.somatic1.variants -o $TUMORID/$TUMORID2.meerkat.somatic2.variants -R $REF/rmsk-hg38.txt -n 1 -D 5 -Q 10 -B $SAMPLEID/$SAMPLEID.recalibrated.bam -I $SAMPLEID/$SAMPLEID.recalibrated.isinfo 
#perl $SCRIPTS/somatic_sv.pl -i $TUMORID/$TUMORID2.meerkat.somatic2.variants -o $TUMORID/$TUMORID2.meerkat.somatic3.variants -R $REF/rmsk-hg38.txt -u 1 -Q 10 -B $SAMPLEID/$SAMPLEID.recalibrated.bam 
#perl $SCRIPTS/somatic_sv.pl -i $TUMORID/$TUMORID2.meerkat.somatic3.variants -o $TUMORID/$TUMORID2.meerkat.somatic4.variants -R $REF/rmsk-hg38.txt -f 1 -Q 10 -B $SAMPLEID/$SAMPLEID.recalibrated.bam 
#perl $SCRIPTS/somatic_sv.pl -i $TUMORID/$TUMORID2.meerkat.somatic4.variants -o $TUMORID/$TUMORID2.meerkat.somatic5.variants -R $REF/rmsk-hg38.txt -e 1 -D 5 -Q 10 -B $TUMORID/$TUMORID.recalibrated.bam -I $TUMORID/$TUMORID.recalibrated.isinfo 
#perl $SCRIPTS/somatic_sv.pl -i $TUMORID/$TUMORID2.meerkat.somatic5.variants -o $TUMORID/$TUMORID2.meerkat.somatic6.variants -R $REF/rmsk-hg38.txt -z 1
#perl $SCRIPTS/somatic_sv.pl -i $TUMORID/$TUMORID2.meerkat.somatic6.variants -o $TUMORID/$TUMORID2.meerkat.somatic.final.variants -R $REF/rmsk-hg38.txt -d 40 -t 20

###Create vcf file - IKKE TILPASSET ENDNU ###
#perl $MER/scripts/meerkat2vcf.pl -i $OUT/somatic_filtered.sample.sorted.nodup.variants -H $MER/Meerkat.example/headerfile -F $REF/Homo_sapiens_assembly38.fasta -o $OUT/somatic_filtered.sample.sorted.nodup.variants.vcf

###Annotation of fusions - IKKE TILPASSET ENDNU ###
#perl $MER/scripts/fusions.pl -i $OUT/somatic_filtered.sample.sorted.nodup.variant -G /data/mark/tools/RAPTR-SV/genes.bed
