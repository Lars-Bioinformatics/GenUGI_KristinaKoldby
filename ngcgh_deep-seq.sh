#!/bin/sh
#
#SBATCH --account sdukoldby_slim      # account
#SBATCH --nodes 1                 # number of nodes
#SBATCH --time 02:00:00            # max time (HH:MM:SS)

WINDOW_SIZE="1000"

# ngCGH -w $WINDOW_SIZE -o ngcgh/ECV2-29-biopsi-C1_tumor_medexome-deep-seq_ngcgh_w$WINDOW_SIZE.txt \
# ECV2-29-blod_normal_tagseq-medexome.connor.recalibrated.bam \
# ECV2-29-biopsi-C1_tumor_tagseq-medexome-deep-seq.connor.recalibrated.bam
#
# ngCGH -w $WINDOW_SIZE -o ngcgh/ECV2-29-biopsi-C2_tumor_medexome-deep-seq_ngcgh_w$WINDOW_SIZE.txt \
# ECV2-29-blod_normal_tagseq-medexome.connor.recalibrated.bam \
# ECV2-29-biopsi-C2_tumor_tagseq-medexome-deep-seq.connor.recalibrated.bam
#
# ngCGH -w $WINDOW_SIZE -o ngcgh/ECV2-29-op-A1_tumor_medexome-deep-seq_ngcgh_w$WINDOW_SIZE.txt \
# ECV2-29-blod_normal_tagseq-medexome.connor.recalibrated.bam \
# ECV2-29-op-A1_tumor_tagseq-medexome-deep-seq.connor.recalibrated.bam
#
# ngCGH -w $WINDOW_SIZE -o ngcgh/ECV2-29-op-A2_tumor_medexome-deep-seq_ngcgh_w$WINDOW_SIZE.txt \
# ECV2-29-blod_normal_tagseq-medexome.connor.recalibrated.bam \
# ECV2-29-op-A2_tumor_tagseq-medexome-deep-seq.connor.recalibrated.bam
#
# ngCGH -w $WINDOW_SIZE -o ngcgh/ECV2-29-op-B1_tumor_medexome-deep-seq_ngcgh_w$WINDOW_SIZE.txt \
# ECV2-29-blod_normal_tagseq-medexome.connor.recalibrated.bam \
# ECV2-29-op-B1_tumor_tagseq-medexome-deep-seq.connor.recalibrated.bam
#
# ngCGH -w $WINDOW_SIZE -o ngcgh/ECV2-29-op-B2_tumor_medexome-deep-seq_ngcgh_w$WINDOW_SIZE.txt \
# ECV2-29-blod_normal_tagseq-medexome.connor.recalibrated.bam \
# ECV2-29-op-B2_tumor_tagseq-medexome-deep-seq.connor.recalibrated.bam
#
# ngCGH -w $WINDOW_SIZE -o ngcgh/ECV2-29-plasma171124_tumor_medexome-deep-seq_ngcgh_w$WINDOW_SIZE.txt \
# ECV2-29-blod_normal_tagseq-medexome.connor.recalibrated.bam \
# ECV2-29-plasma171124_tumor_tagseq-medexome.connor.recalibrated.bam
#
# ngCGH -w $WINDOW_SIZE -o ngcgh/ECV2-29-plasma180119_tumor_medexome-deep-seq_ngcgh_w$WINDOW_SIZE.txt \
# ECV2-29-blod_normal_tagseq-medexome.connor.recalibrated.bam \
# ECV2-29-plasma180119_tumor_tagseq-medexome.connor.recalibrated.bam
#
# ngCGH -w $WINDOW_SIZE -o ngcgh/ECV2-29-plasma180619_tumor_medexome-deep-seq_ngcgh_w$WINDOW_SIZE.txt \
# ECV2-29-blod_normal_tagseq-medexome.connor.recalibrated.bam \
# ECV2-29-plasma180619_tumor_tagseq-medexome.connor.recalibrated.bam
#
# ngCGH -w $WINDOW_SIZE -o ngcgh/ECV2-31-biopsi-F1_tumor_medexome-deep-seq_ngcgh_w$WINDOW_SIZE.txt \
# ECV2-31-blod_normal_tagseq-medexome.connor.recalibrated.bam \
# ECV2-31-biopsi-F1_tumor_tagseq-medexome-deep-seq.connor.recalibrated.bam
#
# ngCGH -w $WINDOW_SIZE -o ngcgh/ECV2-31-biopsi-F2_tumor_medexome-deep-seq_ngcgh_w$WINDOW_SIZE.txt \
# ECV2-31-blod_normal_tagseq-medexome.connor.recalibrated.bam \
# ECV2-31-biopsi-F2_tumor_tagseq-medexome-deep-seq.connor.recalibrated.bam
#
# ngCGH -w $WINDOW_SIZE -o ngcgh/ECV2-31-op-D1_tumor_medexome-deep-seq_ngcgh_w$WINDOW_SIZE.txt \
# ECV2-31-blod_normal_tagseq-medexome.connor.recalibrated.bam \
# ECV2-31-op-D1_tumor_tagseq-medexome-deep-seq.connor.recalibrated.bam
#
# ngCGH -w $WINDOW_SIZE -o ngcgh/ECV2-31-op-D2_tumor_medexome-deep-seq_ngcgh_w$WINDOW_SIZE.txt \
# ECV2-31-blod_normal_tagseq-medexome.connor.recalibrated.bam \
# ECV2-31-op-D2_tumor_tagseq-medexome-deep-seq.connor.recalibrated.bam
#
# ngCGH -w $WINDOW_SIZE -o ngcgh/ECV2-31-op-E1_tumor_medexome-deep-seq_ngcgh_w$WINDOW_SIZE.txt \
# ECV2-31-blod_normal_tagseq-medexome.connor.recalibrated.bam \
# ECV2-31-op-E1_tumor_tagseq-medexome-deep-seq.connor.recalibrated.bam
#
# ngCGH -w $WINDOW_SIZE -o ngcgh/ECV2-31-op-E2_tumor_medexome-deep-seq_ngcgh_w$WINDOW_SIZE.txt \
# ECV2-31-blod_normal_tagseq-medexome.connor.recalibrated.bam \
# ECV2-31-op-E2_tumor_tagseq-medexome-deep-seq.connor.recalibrated.bam
#
# ngCGH -w $WINDOW_SIZE -o ngcgh/ECV2-31-plasma171215_tumor_medexome-deep-seq_ngcgh_w$WINDOW_SIZE.txt \
# ECV2-31-blod_normal_tagseq-medexome.connor.recalibrated.bam \
# ECV2-31-plasma171215_tumor_tagseq-medexome.connor.recalibrated.bam
#
# ngCGH -w $WINDOW_SIZE -o ngcgh/ECV2-31-plasma180202_tumor_medexome-deep-seq_ngcgh_w$WINDOW_SIZE.txt \
# ECV2-31-blod_normal_tagseq-medexome.connor.recalibrated.bam \
# ECV2-31-plasma180202_tumor_tagseq-medexome.connor.recalibrated.bam
#
# ngCGH -w $WINDOW_SIZE -o ngcgh/ECV2-35-biopsi-G1_tumor_medexome-deep-seq_ngcgh_w$WINDOW_SIZE.txt \
# ECV2-35-blod_normal_tagseq-medexome.connor.recalibrated.bam \
# ECV2-35-biopsi-G1_tumor_tagseq-medexome-deep-seq.connor.recalibrated.bam
#
# ngCGH -w $WINDOW_SIZE -o ngcgh/ECV2-35-biopsi-G2_tumor_medexome-deep-seq_ngcgh_w$WINDOW_SIZE.txt \
# ECV2-35-blod_normal_tagseq-medexome.connor.recalibrated.bam \
# ECV2-35-biopsi-G2_tumor_tagseq-medexome-deep-seq.connor.recalibrated.bam
#
# ngCGH -w $WINDOW_SIZE -o ngcgh/ECV2-35-op-01_tumor_medexome-deep-seq_ngcgh_w$WINDOW_SIZE.txt \
# ECV2-35-blod_normal_tagseq-medexome.connor.recalibrated.bam \
# ECV2-35-op-01_tumor_tagseq-medexome-deep-seq.connor.recalibrated.bam
#
# ngCGH -w $WINDOW_SIZE -o ngcgh/ECV2-35-op-02_tumor_medexome-deep-seq_ngcgh_w$WINDOW_SIZE.txt \
# ECV2-35-blod_normal_tagseq-medexome.connor.recalibrated.bam \
# ECV2-35-op-02_tumor_tagseq-medexome-deep-seq.connor.recalibrated.bam
#
# ngCGH -w $WINDOW_SIZE -o ngcgh/ECV2-35-op-03_tumor_medexome-deep-seq_ngcgh_w$WINDOW_SIZE.txt \
# ECV2-35-blod_normal_tagseq-medexome.connor.recalibrated.bam \
# ECV2-35-op-03_tumor_tagseq-medexome-deep-seq.connor.recalibrated.bam
#
# ngCGH -w $WINDOW_SIZE -o ngcgh/ECV2-35-op-04_tumor_medexome-deep-seq_ngcgh_w$WINDOW_SIZE.txt \
# ECV2-35-blod_normal_tagseq-medexome.connor.recalibrated.bam \
# ECV2-35-op-04_tumor_tagseq-medexome-deep-seq.connor.recalibrated.bam
#
# ngCGH -w $WINDOW_SIZE -o ngcgh/ECV2-35-op-06_tumor_medexome-deep-seq_ngcgh_w$WINDOW_SIZE.txt \
# ECV2-35-blod_normal_tagseq-medexome.connor.recalibrated.bam \
# ECV2-35-op-06_tumor_tagseq-medexome-deep-seq.connor.recalibrated.bam
#
# ngCGH -w $WINDOW_SIZE -o ngcgh/ECV2-35-op-07_tumor_medexome-deep-seq_ngcgh_w$WINDOW_SIZE.txt \
# ECV2-35-blod_normal_tagseq-medexome.connor.recalibrated.bam \
# ECV2-35-op-07_tumor_tagseq-medexome-deep-seq.connor.recalibrated.bam
#
# ngCGH -w $WINDOW_SIZE -o ngcgh/ECV2-35-plasma180102_tumor_medexome-deep-seq_ngcgh_w$WINDOW_SIZE.txt \
# ECV2-35-blod_normal_tagseq-medexome.connor.recalibrated.bam \
# ECV2-35-plasma180102_tumor_tagseq-medexome.connor.recalibrated.bam
#
# ngCGH -w $WINDOW_SIZE -o ngcgh/ECV2-35-plasma180316_tumor_medexome-deep-seq_ngcgh_w$WINDOW_SIZE.txt \
# ECV2-35-blod_normal_tagseq-medexome.connor.recalibrated.bam \
# ECV2-35-plasma180316_tumor_tagseq-medexome.connor.recalibrated.bam
#
# ngCGH -w $WINDOW_SIZE -o ngcgh/ECV2-35-plasma180601_tumor_medexome-deep-seq_ngcgh_w$WINDOW_SIZE.txt \
# ECV2-35-blod_normal_tagseq-medexome.connor.recalibrated.bam \
# ECV2-35-plasma180601_tumor_tagseq-medexome.connor.recalibrated.bam

# ngCGH -w $WINDOW_SIZE -o ngcgh/ECV2-4-biopsi-H1_tumor_medexome-deep-seq_ngcgh_w$WINDOW_SIZE.txt \
# ECV2-4-blod_normal_tagseq-medexome.connor.recalibrated.bam \
# ECV2-4-biopsi-H1_tumor_tagseq-medexome-deep-seq.connor.recalibrated.bam
#
# ngCGH -w $WINDOW_SIZE -o ngcgh/ECV2-4-biopsi-H2_tumor_medexome-deep-seq_ngcgh_w$WINDOW_SIZE.txt \
# ECV2-4-blod_normal_tagseq-medexome.connor.recalibrated.bam \
# ECV2-4-biopsi-H2_tumor_tagseq-medexome-deep-seq.connor.recalibrated.bam
#
# ngCGH -w $WINDOW_SIZE -o ngcgh/ECV2-4-op-1401_tumor_medexome-deep-seq_ngcgh_w$WINDOW_SIZE.txt \
# ECV2-4-blod_normal_tagseq-medexome.connor.recalibrated.bam \
# ECV2-4-op-1401_tumor_tagseq-medexome-deep-seq.connor.recalibrated.bam
#
# ngCGH -w $WINDOW_SIZE -o ngcgh/ECV2-4-op-1801_tumor_medexome-deep-seq_ngcgh_w$WINDOW_SIZE.txt \
# ECV2-4-blod_normal_tagseq-medexome.connor.recalibrated.bam \
# ECV2-4-op-1801_tumor_tagseq-medexome-deep-seq.connor.recalibrated.bam
#
# ngCGH -w $WINDOW_SIZE -o ngcgh/ECV2-4-op-1802_tumor_medexome-deep-seq_ngcgh_w$WINDOW_SIZE.txt \
# ECV2-4-blod_normal_tagseq-medexome.connor.recalibrated.bam \
# ECV2-4-op-1802_tumor_tagseq-medexome-deep-seq.connor.recalibrated.bam
#
# ngCGH -w $WINDOW_SIZE -o ngcgh/ECV2-4-plasma170503_tumor_medexome-deep-seq_ngcgh_w$WINDOW_SIZE.txt \
# ECV2-4-blod_normal_tagseq-medexome.connor.recalibrated.bam \
# ECV2-4-plasma170503_tumor_tagseq-medexome.connor.recalibrated.bam
#
# ngCGH -w $WINDOW_SIZE -o ngcgh/ECV2-4-plasma170712_tumor_medexome-deep-seq_ngcgh_w$WINDOW_SIZE.txt \
# ECV2-4-blod_normal_tagseq-medexome.connor.recalibrated.bam \
# ECV2-4-plasma170712_tumor_tagseq-medexome.connor.recalibrated.bam
#
# ngCGH -w $WINDOW_SIZE -o ngcgh/ECV2-4-plasma180205_tumor_medexome-deep-seq_ngcgh_w$WINDOW_SIZE.txt \
# ECV2-4-blod_normal_tagseq-medexome.connor.recalibrated.bam \
# ECV2-4-plasma180205_tumor_tagseq-medexome.connor.recalibrated.bam
#
# ngCGH -w $WINDOW_SIZE -o ngcgh/ECV2-8-biopsi-I1_tumor_medexome-deep-seq_ngcgh_w$WINDOW_SIZE.txt \
# ECV2-8-blod_normal_tagseq-medexome.connor.recalibrated.bam \
# ECV2-8-biopsi-I1_tumor_tagseq-medexome-deep-seq.connor.recalibrated.bam
#
# ngCGH -w $WINDOW_SIZE -o ngcgh/ECV2-8-biopsi-I2_tumor_medexome-deep-seq_ngcgh_w$WINDOW_SIZE.txt \
# ECV2-8-blod_normal_tagseq-medexome.connor.recalibrated.bam \
# ECV2-8-biopsi-I2_tumor_tagseq-medexome-deep-seq.connor.recalibrated.bam
#
# ngCGH -w $WINDOW_SIZE -o ngcgh/ECV2-8-op-01_tumor_medexome-deep-seq_ngcgh_w$WINDOW_SIZE.txt \
# ECV2-8-blod_normal_tagseq-medexome.connor.recalibrated.bam \
# ECV2-8-op-01_tumor_tagseq-medexome-deep-seq.connor.recalibrated.bam
#
# ngCGH -w $WINDOW_SIZE -o ngcgh/ECV2-8-op-02_tumor_medexome-deep-seq_ngcgh_w$WINDOW_SIZE.txt \
# ECV2-8-blod_normal_tagseq-medexome.connor.recalibrated.bam \
# ECV2-8-op-02_tumor_tagseq-medexome-deep-seq.connor.recalibrated.bam
#
# ngCGH -w $WINDOW_SIZE -o ngcgh/ECV2-8-op-03_tumor_medexome-deep-seq_ngcgh_w$WINDOW_SIZE.txt \
# ECV2-8-blod_normal_tagseq-medexome.connor.recalibrated.bam \
# ECV2-8-op-03_tumor_tagseq-medexome-deep-seq.connor.recalibrated.bam
#
# ngCGH -w $WINDOW_SIZE -o ngcgh/ECV2-8-plasma170522_tumor_medexome-deep-seq_ngcgh_w$WINDOW_SIZE.txt \
# ECV2-8-blod_normal_tagseq-medexome.connor.recalibrated.bam \
# ECV2-8-plasma170522_tumor_tagseq-medexome.connor.recalibrated.bam
#
# ngCGH -w $WINDOW_SIZE -o ngcgh/ECV2-8-plasma170728_tumor_medexome-deep-seq_ngcgh_w$WINDOW_SIZE.txt \
# ECV2-8-blod_normal_tagseq-medexome.connor.recalibrated.bam \
# ECV2-8-plasma170728_tumor_tagseq-medexome.connor.recalibrated.bam
#
# ngCGH -w $WINDOW_SIZE -o ngcgh/ECV2-8-plasma180806_tumor_medexome-deep-seq_ngcgh_w$WINDOW_SIZE.txt \
# ECV2-8-blod_normal_tagseq-medexome.connor.recalibrated.bam \
# ECV2-8-plasma180806_tumor_tagseq-medexome.connor.recalibrated.bam
#
# ngCGH -w $WINDOW_SIZE -o ngcgh/PC1-10-op-P1_tumor_medexome-deep-seq_ngcgh_w$WINDOW_SIZE.txt \
# PC1-10-blod_normal_tagseq-medexome.connor.recalibrated.bam \
# PC1-10-op-P1_tumor_tagseq-medexome-deep-seq.connor.recalibrated.bam
#
# ngCGH -w $WINDOW_SIZE -o ngcgh/PC1-10-op-P2_tumor_medexome-deep-seq_ngcgh_w$WINDOW_SIZE.txt \
# PC1-10-blod_normal_tagseq-medexome.connor.recalibrated.bam \
# PC1-10-op-P2_tumor_tagseq-medexome-deep-seq.connor.recalibrated.bam
#
# ngCGH -w $WINDOW_SIZE -o ngcgh/PC1-10-op-P3_tumor_medexome-deep-seq_ngcgh_w$WINDOW_SIZE.txt \
# PC1-10-blod_normal_tagseq-medexome.connor.recalibrated.bam \
# PC1-10-op-P3_tumor_tagseq-medexome-deep-seq.connor.recalibrated.bam
#
# ngCGH -w $WINDOW_SIZE -o ngcgh/PC1-10-op-P4_tumor_medexome-deep-seq_ngcgh_w$WINDOW_SIZE.txt \
# PC1-10-blod_normal_tagseq-medexome.connor.recalibrated.bam \
# PC1-10-op-P4_tumor_tagseq-medexome-deep-seq.connor.recalibrated.bam
#
# ngCGH -w $WINDOW_SIZE -o ngcgh/PC1-10-plasma170928_tumor_medexome-deep-seq_ngcgh_w$WINDOW_SIZE.txt \
# PC1-10-blod_normal_tagseq-medexome.connor.recalibrated.bam \
# PC1-10-plasma170928_tumor_tagseq-medexome.connor.recalibrated.bam
#
# ngCGH -w $WINDOW_SIZE -o ngcgh/PC1-14-op-P5_tumor_medexome-deep-seq_ngcgh_w$WINDOW_SIZE.txt \
# PC1-14-blod_normal_tagseq-medexome.connor.recalibrated.bam \
# PC1-14-op-P5_tumor_tagseq-medexome-deep-seq.connor.recalibrated.bam
#
# ngCGH -w $WINDOW_SIZE -o ngcgh/PC1-14-op-P6_tumor_medexome-deep-seq_ngcgh_w$WINDOW_SIZE.txt \
# PC1-14-blod_normal_tagseq-medexome.connor.recalibrated.bam \
# PC1-14-op-P6_tumor_tagseq-medexome-deep-seq.connor.recalibrated.bam
#
# ngCGH -w $WINDOW_SIZE -o ngcgh/PC1-14-op-P7_tumor_medexome-deep-seq_ngcgh_w$WINDOW_SIZE.txt \
# PC1-14-blod_normal_tagseq-medexome.connor.recalibrated.bam \
# PC1-14-op-P7_tumor_tagseq-medexome-deep-seq.connor.recalibrated.bam
#
# ngCGH -w $WINDOW_SIZE -o ngcgh/PC1-14-op-P8_tumor_medexome-deep-seq_ngcgh_w$WINDOW_SIZE.txt \
# PC1-14-blod_normal_tagseq-medexome.connor.recalibrated.bam \
# PC1-14-op-P8_tumor_tagseq-medexome-deep-seq.connor.recalibrated.bam
#
# ngCGH -w $WINDOW_SIZE -o ngcgh/PC1-14-plasma180111_tumor_medexome-deep-seq_ngcgh_w$WINDOW_SIZE.txt \
# PC1-14-blod_normal_tagseq-medexome.connor.recalibrated.bam \
# PC1-14-plasma180111_tumor_tagseq-medexome.connor.recalibrated.bam
#
# ngCGH -w $WINDOW_SIZE -o ngcgh/PC1-14-plasma180907_tumor_medexome-deep-seq_ngcgh_w$WINDOW_SIZE.txt \
# PC1-14-blod_normal_tagseq-medexome.connor.recalibrated.bam \
# PC1-14-plasma180907_tumor_tagseq-medexome.connor.recalibrated.bam
#
# ngCGH -w $WINDOW_SIZE -o ngcgh/PC1-18-op-P10_tumor_medexome-deep-seq_ngcgh_w$WINDOW_SIZE.txt \
# PC1-18-blod_normal_tagseq-medexome.connor.recalibrated.bam \
# PC1-18-op-P10_tumor_tagseq-medexome-deep-seq.connor.recalibrated.bam
#
# ngCGH -w $WINDOW_SIZE -o ngcgh/PC1-18-op-P11_tumor_medexome-deep-seq_ngcgh_w$WINDOW_SIZE.txt \
# PC1-18-blod_normal_tagseq-medexome.connor.recalibrated.bam \
# PC1-18-op-P11_tumor_tagseq-medexome-deep-seq.connor.recalibrated.bam
#
# ngCGH -w $WINDOW_SIZE -o ngcgh/PC1-18-op-P12_tumor_medexome-deep-seq_ngcgh_w$WINDOW_SIZE.txt \
# PC1-18-blod_normal_tagseq-medexome.connor.recalibrated.bam \
# PC1-18-op-P12_tumor_tagseq-medexome-deep-seq.connor.recalibrated.bam
#
# ngCGH -w $WINDOW_SIZE -o ngcgh/PC1-18-op-P9_tumor_medexome-deep-seq_ngcgh_w$WINDOW_SIZE.txt \
# PC1-18-blod_normal_tagseq-medexome.connor.recalibrated.bam \
# PC1-18-op-P9_tumor_tagseq-medexome-deep-seq.connor.recalibrated.bam
#
# ngCGH -w $WINDOW_SIZE -o ngcgh/PC1-18-plasma180219_tumor_medexome-deep-seq_ngcgh_w$WINDOW_SIZE.txt \
# PC1-18-blod_normal_tagseq-medexome.connor.recalibrated.bam \
# PC1-18-plasma180219_tumor_tagseq-medexome.connor.recalibrated.bam

#######################################
############ MERGED SAMPLES ###########
#######################################

# ngCGH -w $WINDOW_SIZE -o ngcgh/ECV2-29-biopsi-merged_tumor_tagseq-medexome-deep-seq_ngcgh_w$WINDOW_SIZE.txt \
# ECV2-29-blod_normal_tagseq-medexome.connor.recalibrated.bam \
# ECV2-29-biopsi-merged_tumor_tagseq-medexome-deep-seq.connor.recalibrated.bam
#
# ngCGH -w $WINDOW_SIZE -o ngcgh/ECV2-29-opA-merged_tumor_tagseq-medexome-deep-seq_ngcgh_w$WINDOW_SIZE.txt \
# ECV2-29-blod_normal_tagseq-medexome.connor.recalibrated.bam \
# ECV2-29-opA-merged_tumor_tagseq-medexome-deep-seq.connor.recalibrated.bam
#
# ngCGH -w $WINDOW_SIZE -o ngcgh/ECV2-29-opB-merged_tumor_tagseq-medexome-deep-seq_ngcgh_w$WINDOW_SIZE.txt \
# ECV2-29-blod_normal_tagseq-medexome.connor.recalibrated.bam \
# ECV2-29-opB-merged_tumor_tagseq-medexome-deep-seq.connor.recalibrated.bam
#
# ngCGH -w $WINDOW_SIZE -o ngcgh/ECV2-31-biopsi-merged_tumor_tagseq-medexome-deep-seq_ngcgh_w$WINDOW_SIZE.txt \
# ECV2-31-blod_normal_tagseq-medexome.connor.recalibrated.bam \
# ECV2-31-biopsi-merged_tumor_tagseq-medexome-deep-seq.connor.recalibrated.bam
#
# ngCGH -w $WINDOW_SIZE -o ngcgh/ECV2-31-opD-merged_tumor_tagseq-medexome-deep-seq_ngcgh_w$WINDOW_SIZE.txt \
# ECV2-31-blod_normal_tagseq-medexome.connor.recalibrated.bam \
# ECV2-31-opD-merged_tumor_tagseq-medexome-deep-seq.connor.recalibrated.bam
#
# ngCGH -w $WINDOW_SIZE -o ngcgh/ECV2-31-opE-merged_tumor_tagseq-medexome-deep-seq_ngcgh_w$WINDOW_SIZE.txt \
# ECV2-31-blod_normal_tagseq-medexome.connor.recalibrated.bam \
# ECV2-31-opE-merged_tumor_tagseq-medexome-deep-seq.connor.recalibrated.bam
#
# ngCGH -w $WINDOW_SIZE -o ngcgh/ECV2-35-biopsi-merged_tumor_tagseq-medexome-deep-seq_ngcgh_w$WINDOW_SIZE.txt \
# ECV2-35-blod_normal_tagseq-medexome.connor.recalibrated.bam \
# ECV2-35-biopsi-merged_tumor_tagseq-medexome-deep-seq.connor.recalibrated.bam
#
# ngCGH -w $WINDOW_SIZE -o ngcgh/ECV2-35-op1to4-merged_tumor_tagseq-medexome-deep-seq_ngcgh_w$WINDOW_SIZE.txt \
# ECV2-35-blod_normal_tagseq-medexome.connor.recalibrated.bam \
# ECV2-35-op1to4-merged_tumor_tagseq-medexome-deep-seq.connor.recalibrated.bam
#
# ngCGH -w $WINDOW_SIZE -o ngcgh/ECV2-35-op6to7-merged_tumor_tagseq-medexome-deep-seq_ngcgh_w$WINDOW_SIZE.txt \
# ECV2-35-blod_normal_tagseq-medexome.connor.recalibrated.bam \
# ECV2-35-op6to7-merged_tumor_tagseq-medexome-deep-seq.connor.recalibrated.bam
#
# ngCGH -w $WINDOW_SIZE -o ngcgh/ECV2-4-biopsi-merged_tumor_tagseq-medexome-deep-seq_ngcgh_w$WINDOW_SIZE.txt \
# ECV2-4-blod_normal_tagseq-medexome.connor.recalibrated.bam \
# ECV2-4-biopsi-merged_tumor_tagseq-medexome-deep-seq.connor.recalibrated.bam
#
# ngCGH -w $WINDOW_SIZE -o ngcgh/ECV2-4-op18-merged_tumor_tagseq-medexome-deep-seq_ngcgh_w$WINDOW_SIZE.txt \
# ECV2-4-blod_normal_tagseq-medexome.connor.recalibrated.bam \
# ECV2-4-op18-merged_tumor_tagseq-medexome-deep-seq.connor.recalibrated.bam
#
# ngCGH -w $WINDOW_SIZE -o ngcgh/ECV2-4-op-merged_tumor_tagseq-medexome-deep-seq_ngcgh_w$WINDOW_SIZE.txt \
# ECV2-4-blod_normal_tagseq-medexome.connor.recalibrated.bam \
# ECV2-4-op-merged_tumor_tagseq-medexome-deep-seq.connor.recalibrated.bam
#
# ngCGH -w $WINDOW_SIZE -o ngcgh/ECV2-8-biopsi-merged_tumor_tagseq-medexome-deep-seq_ngcgh_w$WINDOW_SIZE.txt \
# ECV2-8-blod_normal_tagseq-medexome.connor.recalibrated.bam \
# ECV2-8-biopsi-merged_tumor_tagseq-medexome-deep-seq.connor.recalibrated.bam
#
# ngCGH -w $WINDOW_SIZE -o ngcgh/ECV2-8-op-merged_tumor_tagseq-medexome-deep-seq_ngcgh_w$WINDOW_SIZE.txt \
# ECV2-8-blod_normal_tagseq-medexome.connor.recalibrated.bam \
# ECV2-8-op-merged_tumor_tagseq-medexome-deep-seq.connor.recalibrated.bam
#
# ngCGH -w $WINDOW_SIZE -o ngcgh/PC1-10-op-merged_tumor_tagseq-medexome-deep-seq_ngcgh_w$WINDOW_SIZE.txt \
# PC1-10-blod_normal_tagseq-medexome.connor.recalibrated.bam \
# PC1-10-op-merged_tumor_tagseq-medexome-deep-seq.connor.recalibrated.bam
#
# ngCGH -w $WINDOW_SIZE -o ngcgh/PC1-14-op-merged_tumor_tagseq-medexome-deep-seq_ngcgh_w$WINDOW_SIZE.txt \
# PC1-14-blod_normal_tagseq-medexome.connor.recalibrated.bam \
# PC1-14-op-merged_tumor_tagseq-medexome-deep-seq.connor.recalibrated.bam
#
ngCGH -w $WINDOW_SIZE -o ngcgh/PC1-18-op-merged_tumor_tagseq-medexome-deep-seq_ngcgh_w$WINDOW_SIZE.txt \
PC1-18-blod_normal_tagseq-medexome.connor.recalibrated.bam \
PC1-18-op-merged_tumor_tagseq-medexome-deep-seq.connor.recalibrated.bam
