#!/bin/sh
#
#SBATCH --account sduhumac_slim      # account
#SBATCH --nodes 1                 # number of nodes
#SBATCH --time 24:00:00            # max time (HH:MM:SS)
#SBATCH --output=/work/sduhumac/kristina/scripts/slurm_log/map_sort_markduplicates_%j.out

# Call script by sbatch map_sort_markduplicates.sh SAMPLEID

set -x
$PWD

### Specify paths to requiered software, input and output folders
SAMPLEID=$1
DAT="/gpfs/gss1/work/sduhumac/kristina/data/genova/plasma"
OUTPUT="/gpfs/gss1/scratch/sduhumac/kristina/data/genova/plasma/output"
REF="/gpfs/gss1/work/sduhumac/tools/ref"
BBMAP="/gpfs/gss1/work/sduhumac/kristina/tools/miniconda3/pkgs/bbmap-38.22-h470a237_0/"
BREAKDANCER="/gpfs/gss1/work/sduhumac/kristina/tools/miniconda3/pkgs/breakdancer-1.4.5-2/bin"


###### Quality filtering
$BBMAP/bin/bbduk.sh -Xmx20g \
in1=/$DAT/${SAMPLEID}_R1.fastq.gz \
in2=/$DAT/${SAMPLEID}_R2.fastq.gz \
out1=/$DAT/${SAMPLEID}_clean_R1.fastq.gz \
out2=/$DAT/${SAMPLEID}_clean_R2.fastq.gz \
ref=$BBMAP/opt/bbmap-38.22-0/resources/adapters.fa \
ktrim=r \
ktrim=l \
k=23 \
mink=11 \
hdist=1 tpe tbo \
qtrim="rl" \
trimq=10 \
maq=10 \
minlen=25

###### BWA alignment
bwa mem -t 24 -M -R '@RG\tID:${SAMPLEID}.lib.run\tLB:${SAMPLEID}.lib\tPL:ILLUMINA\tSM:${SAMPLEID}' $REF/hg19.fa $DAT/${SAMPLEID}_clean_R1.fastq.gz $DAT/${SAMPLEID}_clean_R2.fastq.gz | samtools sort - > $OUTPUT/${SAMPLEID}_sorted.bam 

### Samtools processing of aligned reads
samtools index $OUTPUT/${SAMPLEID}_sorted.bam $OUTPUT/${SAMPLEID}_sorted.bai

### PCR duplicate removal
#picard MarkDuplicates \
#INPUT=$OUTPUT/${SAMPLEID}_sorted.bam \
#OUTPUT=$OUTPUT/${SAMPLEID}_sorted_nodup.bam \
#METRICS_FILE=${SAMPLEID}_dup.metrics \
#REMOVE_DUPLICATES=TRUE \
#VALIDATION_STRINGENCY=LENIENT \
#CREATE_INDEX=true

### Remove old files
#rm -f $OUTPUT/${SAMPLEID}_sorted.bam
#rm -f $OUTPUT/${SAMPLEID}_sorted.bai

### Quality control with qualimap
#qualimap bamqc -bam $OUTPUT/${SAMPLEID}_sorted_nodup.bam -c -outdir $OUTPUT/qc/$SAMPLEID -nt 24 --java-mem-size=20G

