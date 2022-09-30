#!/bin/sh
#
#SBATCH --account sduhumac_fat      # account
#SBATCH --nodes 1                 # number of nodes
#SBATCH --time 24:00:00            # max time (HH:MM:SS)

### Specify paths to requiered software, input and output folders
DAT="/gpfs/gss1/work/sduhumac/kristina/data/genova"

###### Compressing unpacked fastq files
cd $DAT
for i in *_R2.fastq;
do
newfile=$(basename $i _R2.fastq)
gzip $DAT/${newfile}_R1.fastq
gzip $DAT/${newfile}_R2.fastq
done