#!/bin/sh
#
#SBATCH --account sduhumac_fat      # account
#SBATCH --nodes 1                 # number of nodes
#SBATCH --time 24:00:00            # max time (HH:MM:SS)
#SBATCH --output=/work/sduhumac/kristina/scripts/slurm_log/merge_bam_markduplicates_%j.out

set -x
$PWD


### Specify paths to requiered software, input and output folders
SAMPLEID=$1
SAMPLEID_1=$2
SAMPLEID_2=$3
OUTPUT="/scratch/sduhumac/kristina/data/genova/plasma/output"
LOGFILE="/scratch/sduhumac/kristina/data/genova/plasma/output/logfiles/$SAMPLEID_log_merge_bam_markduplicates.txt" 

### Merge bam files
if ! [[ $2 == '' ]]; then
  samtools merge $OUTPUT/${SAMPLEID}_sorted.bam $OUTPUT/${SAMPLEID_1}_sorted.bam $OUTPUT/${SAMPLEID_2}_sorted.bam
fi

### PCR duplicate removal
### Kan kræve mange RAM
picard MarkDuplicates -Xmx40g \
INPUT=$OUTPUT/${SAMPLEID}_sorted.bam \
OUTPUT=$OUTPUT/${SAMPLEID}_sorted_nodup.bam \
METRICS_FILE=$OUTPUT/${SAMPLEID}_dup.metrics \
REMOVE_DUPLICATES=TRUE \
VALIDATION_STRINGENCY=LENIENT \
CREATE_INDEX=true


### Remove old files
#rm -f $OUTPUT/${SAMPLEID}_sorted.bam
#rm -f $OUTPUT/${SAMPLEID}_sorted.bai

### Quality control with qualimap
qualimap bamqc \
-bam $OUTPUT/${SAMPLEID}_sorted_nodup.bam \
-nt 24 \
-c \
-sd \
-outdir $OUTPUT/qc/$SAMPLEID/ \
--java-mem-size=20G

### Quality control with fastqc
fastqc $OUTPUT/${SAMPLEID}_sorted_nodup.bam -o $OUTPUT/qc/fastqc/
