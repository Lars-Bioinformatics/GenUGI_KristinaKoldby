#!/bin/sh
#
#SBATCH --account sdukoldby_fat      # account
#SBATCH --nodes 1                 # number of nodes
#SBATCH --time 24:00:00            # max time (HH:MM:SS)

### RUN FROM FOLDER WITH BAM-FILES ###

REF="/work/sdukoldby/resources/hg38/Homo_sapiens_assembly38.fasta"
EXCLUDE="/work/sdukoldby/resources/hg38/exclude.cnvnator_100bp.GRCh38.20170403.bed"
OUTDIR="lumpy"
OUTDIR_JOINT="lumpy/joint_calling_all"
OUTDIR_NORMALS="lumpy/joint_calling_normals"

# smoove call -x --name G45-ECV2-29-biopsi-C1 --outdir $OUTDIR --exclude $EXCLUDE --fasta $REF -p 24 --genotype G45-ECV2-29-biopsi-C1_truseq-nano-genome_HGL2LDSXX_S21.recalibrated.bam G45-ECV2-29-blod_truseq-nano-genome_HGL2LDSXX_S30.recalibrated.bam
# smoove call -x --name G45-ECV2-31-biopsi-F1 --outdir $OUTDIR --exclude $EXCLUDE --fasta $REF -p 24 --genotype G45-ECV2-31-biopsi-F1_truseq-nano-genome_HGL2LDSXX_S22.recalibrated.bam G45-ECV2-31-blod_truseq-nano-genome_HGL2LDSXX_S27.recalibrated.bam
# smoove call -x --name G45-ECV2-35-biopsi-G1 --outdir $OUTDIR --exclude $EXCLUDE --fasta $REF -p 24 --genotype G45-ECV2-35-biopsi-G1_truseq-nano-genome_HGL2LDSXX_S23.recalibrated.bam G45-ECV2-35-blod_truseq-nano-genome_HGL2LDSXX_S28.recalibrated.bam
# smoove call -x --name G45-ECV2-4-biopsi-H2 --outdir $OUTDIR --exclude $EXCLUDE --fasta $REF -p 24 --genotype G45-ECV2-4-biopsi-H2_truseq-nano-genome_HGL2LDSXX_S24.recalibrated.bam G45-ECV2-4-blod_truseq-nano-genome_HGL2LDSXX_S26.recalibrated.bam
# smoove call -x --name G45-ECV2-8-biopsi-I2 --outdir $OUTDIR --exclude $EXCLUDE --fasta $REF -p 24 --genotype G45-ECV2-8-biopsi-I2_truseq-nano-genome_HGL2LDSXX_S25.recalibrated.bam G45-ECV2-8-blod_truseq-nano-genome_HGL2LDSXX_S29.recalibrated.bam

### Joint calling all samples ###
smoove call -x --name G45-ECV2-joint-all --outdir $OUTDIR_JOINT --exclude $EXCLUDE --fasta $REF -p 24 \
--genotype G45-ECV2-29-biopsi-C1_truseq-nano-genome_HGL2LDSXX_S21.recalibrated.bam \
G45-ECV2-29-blod_truseq-nano-genome_HGL2LDSXX_S30.recalibrated.bam \
G45-ECV2-31-biopsi-F1_truseq-nano-genome_HGL2LDSXX_S22.recalibrated.bam \
G45-ECV2-31-blod_truseq-nano-genome_HGL2LDSXX_S27.recalibrated.bam \
G45-ECV2-35-biopsi-G1_truseq-nano-genome_HGL2LDSXX_S23.recalibrated.bam \
G45-ECV2-35-blod_truseq-nano-genome_HGL2LDSXX_S28.recalibrated.bam \
G45-ECV2-4-biopsi-H2_truseq-nano-genome_HGL2LDSXX_S24.recalibrated.bam \
G45-ECV2-4-blod_truseq-nano-genome_HGL2LDSXX_S26.recalibrated.bam \
G45-ECV2-8-biopsi-I2_truseq-nano-genome_HGL2LDSXX_S25.recalibrated.bam \
G45-ECV2-8-blod_truseq-nano-genome_HGL2LDSXX_S29.recalibrated.bam


### Joint calling all normals ###
# smoove call -x --name G45-ECV2-joint --outdir $OUTDIR_NORMALS --exclude $EXCLUDE --fasta $REF -p 24 \
# --genotype G45-ECV2-29-blod_truseq-nano-genome_HGL2LDSXX_S30.recalibrated.bam \
# G45-ECV2-31-blod_truseq-nano-genome_HGL2LDSXX_S27.recalibrated.bam \
# G45-ECV2-35-blod_truseq-nano-genome_HGL2LDSXX_S28.recalibrated.bam \
# G45-ECV2-4-blod_truseq-nano-genome_HGL2LDSXX_S26.recalibrated.bam \
# G45-ECV2-8-blod_truseq-nano-genome_HGL2LDSXX_S29.recalibrated.bam
