#!/bin/sh
#
#SBATCH --account sduhumac_fat   # account
#SBATCH --nodes 1                 # number of nodes
#SBATCH --time 24:00:00           # max time (HH:MM:SS)

### Produce a GC Wiggle track file ###
#sequenza-utils gc_wiggle -w 50 -f Homo_sapiens_assembly38.fasta -o Homo_sapiens_assembly38.gc50Base.wig.gz

### Process BAM and Wiggle files to produce a seqz file ###
# sequenza-utils bam2seqz \
#   -n ECV2-4-blod_normal_tagseq-medexome.connor.recalibrated.bam \
#   -t ECV2-4-biopsi-H1_tumor_tagseq-medexome.connor.recalibrated.bam \
#   --fasta Homo_sapiens_assembly38.fasta \
#   -gc Homo_sapiens_assembly38.gc50Base.wig.gz \
#   -o ECV2-4-biopsi-H1.seqz.gz

### Post-process by binning the original seqz file ###
sequenza-utils seqz_binning --seqz ECV2-4-biopsi-H1.seqz.gz -w 50 -o ECV2-4-biopsi-H1.bins.seqz.gz
