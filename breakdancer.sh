#!/bin/sh
#
#SBATCH --account sduhumac_fat      # account
#SBATCH --nodes 1                 # number of nodes
#SBATCH --time 24:00:00            # max time (HH:MM:SS)


SAMPLEID=$1
INPUT="/work/sdukoldby/data/G45-2016_genugi/190325_A00653_0016_AHGL2LDSXX/BaseCalls/hg38"
OUTPUT="/work/sdukoldby/data/G45-2016_genugi/190325_A00653_0016_AHGL2LDSXX/BaseCalls/hg38/breakdancer"
SRC="/work/sdukoldby/tools/miniconda3/pkgs/breakdancer-1.4.5-2/bin"

### Run breakdancer
#cpanm --local-lib=~/perl5 local::lib && eval $(perl -I ~/perl5/lib/perl5/ -Mlocal::lib)
$SRC/bam2cfg.pl $INPUT/${SAMPLEID}.recalibrated.bam  > $OUTPUT/${SAMPLEID}.cfg
breakdancer-max $OUTPUT/${SAMPLEID}.cfg -d -g -h > $OUTPUT/${SAMPLEID}.ctx

