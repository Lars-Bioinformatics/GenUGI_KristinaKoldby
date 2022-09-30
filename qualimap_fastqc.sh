#!/bin/sh
#
#SBATCH --account sduhumac_slim      # account
#SBATCH --nodes 1                 # number of nodes
#SBATCH --time 24:00:00            # max time (HH:MM:SS)
#SBATCH --output=/work/sduhumac/kristina/scripts/slurm_log/qualimap_fastqc_%j.out

set -x
echo $PWD


### Specify paths to requiered software, input and output folders
SAMPLEID=$1
OUTPUT="/scratch/sduhumac/kristina/data/genova/plasma/output"
#OUTPUT="/work/sduhumac/kristina/data/genova/tumor_blood/output"

### Quality control with qualimap
qualimap bamqc \
-bam $OUTPUT/${SAMPLEID}_sorted_nodup.bam \
-nt 24 \
-c \
-sd \
-outdir $OUTPUT/qc/$SAMPLEID/ \
--java-mem-size=20G

### Quality control with fastqc
#fastqc $OUTPUT/${SAMPLEID}_sorted_nodup.bam -o $OUTPUT/qc/fastqc/
