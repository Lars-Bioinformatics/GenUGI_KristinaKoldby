#!/bin/sh
#
#SBATCH --account sdukoldby_fat
#SBATCH --nodes 1
#SBATCH --time 05:00:00

# cp ECV2-29-blod_normal_takara-tagseq_BHGTHHDSXX_R1_001.fastq.gz 	   exome_fastq_merged/ECV2-29-blod_normal_tagseq-medexome_R1_001.fastq.gz &
# cp ECV2-29-plasma171124_tumor_takara-tagseq_BHGTHHDSXX_R1_001.fastq.gz exome_fastq_merged/ECV2-29-plasma171124_tumor_tagseq-medexome_R1_001.fastq.gz &
# cp ECV2-29-plasma180119_tumor_takara-tagseq_BHGTHHDSXX_R1_001.fastq.gz exome_fastq_merged/ECV2-29-plasma180119_tumor_tagseq-medexome_R1_001.fastq.gz &
# cp ECV2-29-plasma180619_tumor_takara-tagseq_BHGTHHDSXX_R1_001.fastq.gz exome_fastq_merged/ECV2-29-plasma180619_tumor_tagseq-medexome_R1_001.fastq.gz &
# cp ECV2-31-biopsi-F1_tumor_takara-tagseq_BHGTHHDSXX_R1_001.fastq.gz    exome_fastq_merged/ECV2-31-biopsi-F1_tumor_tagseq-medexome_R1_001.fastq.gz &
# cp ECV2-31-biopsi-F2_tumor_takara-tagseq_BHGTHHDSXX_R1_001.fastq.gz    exome_fastq_merged/ECV2-31-biopsi-F2_tumor_tagseq-medexome_R1_001.fastq.gz &
# cp ECV2-31-blod_normal_takara-tagseq_BHGTHHDSXX_R1_001.fastq.gz        exome_fastq_merged/ECV2-31-blod_normal_tagseq-medexome_R1_001.fastq.gz &
# cp ECV2-31-op-D1_tumor_takara-tagseq_BHGTHHDSXX_R1_001.fastq.gz        exome_fastq_merged/ECV2-31-op-D1_tumor_tagseq-medexome_R1_001.fastq.gz &
# cp ECV2-31-op-E1_tumor_takara-tagseq_BHGTHHDSXX_R1_001.fastq.gz        exome_fastq_merged/ECV2-31-op-E1_tumor_tagseq-medexome_R1_001.fastq.gz &
# cp ECV2-31-op-E2_tumor_takara-tagseq_BHGTHHDSXX_R1_001.fastq.gz        exome_fastq_merged/ECV2-31-op-E2_tumor_tagseq-medexome_R1_001.fastq.gz &
# cp ECV2-31-plasma171215_tumor_takara-tagseq_BHGTHHDSXX_R1_001.fastq.gz exome_fastq_merged/ECV2-31-plasma171215_tumor_tagseq-medexome_R1_001.fastq.gz &
# cp ECV2-31-plasma180202_tumor_takara-tagseq_BHGTHHDSXX_R1_001.fastq.gz exome_fastq_merged/ECV2-31-plasma180202_tumor_tagseq-medexome_R1_001.fastq.gz &
# cp ECV2-35-biopsi-G1_tumor_takara-tagseq_BHGTHHDSXX_R1_001.fastq.gz    exome_fastq_merged/ECV2-35-biopsi-G1_tumor_tagseq-medexome_R1_001.fastq.gz &
# cp ECV2-35-biopsi-G2_tumor_takara-tagseq_BHGTHHDSXX_R1_001.fastq.gz    exome_fastq_merged/ECV2-35-biopsi-G2_tumor_tagseq-medexome_R1_001.fastq.gz &
# cp ECV2-35-plasma180102_tumor_takara-tagseq_BHGTHHDSXX_R1_001.fastq.gz exome_fastq_merged/ECV2-35-plasma180102_tumor_tagseq-medexome_R1_001.fastq.gz &
# cp ECV2-35-plasma180316_tumor_takara-tagseq_BHGTHHDSXX_R1_001.fastq.gz exome_fastq_merged/ECV2-35-plasma180316_tumor_tagseq-medexome_R1_001.fastq.gz &
# cp ECV2-35-plasma180601_tumor_takara-tagseq_BHGTHHDSXX_R1_001.fastq.gz exome_fastq_merged/ECV2-35-plasma180601_tumor_tagseq-medexome_R1_001.fastq.gz &
# cp ECV2-35-recidiv_takara-tagseq_BHGTHHDSXX_R1_001.fastq.gz            exome_fastq_merged/ECV2-35-recidiv_tagseq-medexome_R1_001.fastq.gz &
# cp ECV2-4-biopsi-H1_tumor_takara-tagseq_BHGTHHDSXX_R1_001.fastq.gz     exome_fastq_merged/ECV2-4-biopsi-H1_tumor_tagseq-medexome_R1_001.fastq.gz &
# cp ECV2-4-biopsi-H2_tumor_takara-tagseq_BHGTHHDSXX_R1_001.fastq.gz     exome_fastq_merged/ECV2-4-biopsi-H2_tumor_tagseq-medexome_R1_001.fastq.gz &
# cp ECV2-4-blod_normal_takara-tagseq_BHGTHHDSXX_R1_001.fastq.gz         exome_fastq_merged/ECV2-4-blod_normal_tagseq-medexome_R1_001.fastq.gz &
# cp ECV2-4-op-1401_tumor_takara-tagseq_BHGTHHDSXX_R1_001.fastq.gz       exome_fastq_merged/ECV2-4-op-1401_tumor_tagseq-medexome_R1_001.fastq.gz &
# cp ECV2-4-op-1801_tumor_takara-tagseq_BHGTHHDSXX_R1_001.fastq.gz       exome_fastq_merged/ECV2-4-op-1801_tumor_tagseq-medexome_R1_001.fastq.gz &
# cp ECV2-4-plasma170712_tumor_takara-tagseq_BHGTHHDSXX_R1_001.fastq.gz  exome_fastq_merged/ECV2-4-plasma170712_tumor_tagseq-medexome_R1_001.fastq.gz &
# cp ECV2-4-plasma180205_tumor_takara-tagseq_BHGTHHDSXX_R1_001.fastq.gz  exome_fastq_merged/ECV2-4-plasma180205_tumor_tagseq-medexome_R1_001.fastq.gz &
# cp ECV2-8-biopsi-I1_tumor_takara-tagseq_BHGTHHDSXX_R1_001.fastq.gz     exome_fastq_merged/ECV2-8-biopsi-I1_tumor_tagseq-medexome_R1_001.fastq.gz &
# cp ECV2-8-biopsi-I2_tumor_takara-tagseq_BHGTHHDSXX_R1_001.fastq.gz     exome_fastq_merged/ECV2-8-biopsi-I2_tumor_tagseq-medexome_R1_001.fastq.gz &
# cp ECV2-8-blod_normal_takara-tagseq_BHGTHHDSXX_R1_001.fastq.gz         exome_fastq_merged/ECV2-8-blod_normal_tagseq-medexome_R1_001.fastq.gz &
# cp ECV2-8-op-01_tumor_takara-tagseq_BHGTHHDSXX_R1_001.fastq.gz         exome_fastq_merged/ECV2-8-op-01_tumor_tagseq-medexome_R1_001.fastq.gz &
# cp ECV2-8-op-02_tumor_takara-tagseq_BHGTHHDSXX_R1_001.fastq.gz         exome_fastq_merged/ECV2-8-op-02_tumor_tagseq-medexome_R1_001.fastq.gz &
# cp ECV2-8-plasma170522_tumor_takara-tagseq_BHGTHHDSXX_R1_001.fastq.gz  exome_fastq_merged/ECV2-8-plasma170522_tumor_tagseq-medexome_R1_001.fastq.gz &
# cp ECV2-8-plasma170728_tumor_takara-tagseq_BHGTHHDSXX_R1_001.fastq.gz  exome_fastq_merged/ECV2-8-plasma170728_tumor_tagseq-medexome_R1_001.fastq.gz &
# cp ECV2-8-plasma180806_tumor_takara-tagseq_BHGTHHDSXX_R1_001.fastq.gz  exome_fastq_merged/ECV2-8-plasma180806_tumor_tagseq-medexome_R1_001.fastq.gz &
# cp PC1-10-blod_normal_takara-tagseq_BHGTHHDSXX_R1_001.fastq.gz         exome_fastq_merged/PC1-10-blod_normal_tagseq-medexome_R1_001.fastq.gz &
# cp PC1-10-op-P1_tumor_takara-tagseq_BHGTHHDSXX_R1_001.fastq.gz         exome_fastq_merged/PC1-10-op-P1_tumor_tagseq-medexome_R1_001.fastq.gz &
# cp PC1-10-op-P2_tumor_takara-tagseq_BHGTHHDSXX_R1_001.fastq.gz         exome_fastq_merged/PC1-10-op-P2_tumor_tagseq-medexome_R1_001.fastq.gz &
# cp PC1-10-op-P3_tumor_takara-tagseq_BHGTHHDSXX_R1_001.fastq.gz         exome_fastq_merged/PC1-10-op-P3_tumor_tagseq-medexome_R1_001.fastq.gz &
# cp PC1-10-op-P4_tumor_takara-tagseq_BHGTHHDSXX_R1_001.fastq.gz         exome_fastq_merged/PC1-10-op-P4_tumor_tagseq-medexome_R1_001.fastq.gz &
# cp PC1-10-plasma170928_tumor_takara-tagseq_BHGTHHDSXX_R1_001.fastq.gz  exome_fastq_merged/PC1-10-plasma170928_tumor_tagseq-medexome_R1_001.fastq.gz &
# cp PC1-14-blod_normal_takara-tagseq_BHGTHHDSXX_R1_001.fastq.gz         exome_fastq_merged/PC1-14-blod_normal_tagseq-medexome_R1_001.fastq.gz &
# cp PC1-14-op-P5_tumor_takara-tagseq_BHGTHHDSXX_R1_001.fastq.gz         exome_fastq_merged/PC1-14-op-P5_tumor_tagseq-medexome_R1_001.fastq.gz &
# cp PC1-14-op-P6_tumor_takara-tagseq_BHGTHHDSXX_R1_001.fastq.gz         exome_fastq_merged/PC1-14-op-P6_tumor_tagseq-medexome_R1_001.fastq.gz &
# cp PC1-14-op-P7_tumor_takara-tagseq_BHGTHHDSXX_R1_001.fastq.gz         exome_fastq_merged/PC1-14-op-P7_tumor_tagseq-medexome_R1_001.fastq.gz &
# cp PC1-14-op-P8_tumor_takara-tagseq_BHGTHHDSXX_R1_001.fastq.gz         exome_fastq_merged/PC1-14-op-P8_tumor_tagseq-medexome_R1_001.fastq.gz &
# cp PC1-14-plasma180111_tumor_takara-tagseq_BHGTHHDSXX_R1_001.fastq.gz  exome_fastq_merged/PC1-14-plasma180111_tumor_tagseq-medexome_R1_001.fastq.gz &
# cp PC1-14-plasma180907_tumor_takara-tagseq_BHGTHHDSXX_R1_001.fastq.gz  exome_fastq_merged/PC1-14-plasma180907_tumor_tagseq-medexome_R1_001.fastq.gz &
# cp PC1-14-recidiv_takara-tagseq_BHGTHHDSXX_R1_001.fastq.gz             exome_fastq_merged/PC1-14-recidiv_tagseq-medexome_R1_001.fastq.gz &
# cp PC1-18-blod_normal_takara-tagseq_BHGTHHDSXX_R1_001.fastq.gz         exome_fastq_merged/PC1-18-blod_normal_tagseq-medexome_R1_001.fastq.gz &
# cp PC1-18-op-P10_tumor_takara-tagseq_BHGTHHDSXX_R1_001.fastq.gz        exome_fastq_merged/PC1-18-op-P10_tumor_tagseq-medexome_R1_001.fastq.gz &
# cp PC1-18-op-P11_tumor_takara-tagseq_BHGTHHDSXX_R1_001.fastq.gz        exome_fastq_merged/PC1-18-op-P11_tumor_tagseq-medexome_R1_001.fastq.gz &
# cp PC1-18-op-P9_tumor_takara-tagseq_BHGTHHDSXX_R1_001.fastq.gz         exome_fastq_merged/PC1-18-op-P9_tumor_tagseq-medexome_R1_001.fastq.gz &
# cp PC1-18-plasma180219_tumor_takara-tagseq_BHGTHHDSXX_R1_001.fastq.gz  exome_fastq_merged/PC1-18-plasma180219_tumor_tagseq-medexome_R1_001.fastq.gz &
#
#
# cp ECV2-29-blod_normal_takara-tagseq_BHGTHHDSXX_R2_001.fastq.gz 	   exome_fastq_merged/ECV2-29-blod_normal_tagseq-medexome_R2_001.fastq.gz &
# cp ECV2-29-plasma171124_tumor_takara-tagseq_BHGTHHDSXX_R2_001.fastq.gz exome_fastq_merged/ECV2-29-plasma171124_tumor_tagseq-medexome_R2_001.fastq.gz &
# cp ECV2-29-plasma180119_tumor_takara-tagseq_BHGTHHDSXX_R2_001.fastq.gz exome_fastq_merged/ECV2-29-plasma180119_tumor_tagseq-medexome_R2_001.fastq.gz &
# cp ECV2-29-plasma180619_tumor_takara-tagseq_BHGTHHDSXX_R2_001.fastq.gz exome_fastq_merged/ECV2-29-plasma180619_tumor_tagseq-medexome_R2_001.fastq.gz &
# cp ECV2-31-biopsi-F1_tumor_takara-tagseq_BHGTHHDSXX_R2_001.fastq.gz    exome_fastq_merged/ECV2-31-biopsi-F1_tumor_tagseq-medexome_R2_001.fastq.gz &
# cp ECV2-31-biopsi-F2_tumor_takara-tagseq_BHGTHHDSXX_R2_001.fastq.gz    exome_fastq_merged/ECV2-31-biopsi-F2_tumor_tagseq-medexome_R2_001.fastq.gz &
# cp ECV2-31-blod_normal_takara-tagseq_BHGTHHDSXX_R2_001.fastq.gz        exome_fastq_merged/ECV2-31-blod_normal_tagseq-medexome_R2_001.fastq.gz &
# cp ECV2-31-op-D1_tumor_takara-tagseq_BHGTHHDSXX_R2_001.fastq.gz        exome_fastq_merged/ECV2-31-op-D1_tumor_tagseq-medexome_R2_001.fastq.gz &
# cp ECV2-31-op-E1_tumor_takara-tagseq_BHGTHHDSXX_R2_001.fastq.gz        exome_fastq_merged/ECV2-31-op-E1_tumor_tagseq-medexome_R2_001.fastq.gz &
# cp ECV2-31-op-E2_tumor_takara-tagseq_BHGTHHDSXX_R2_001.fastq.gz        exome_fastq_merged/ECV2-31-op-E2_tumor_tagseq-medexome_R2_001.fastq.gz &
# cp ECV2-31-plasma171215_tumor_takara-tagseq_BHGTHHDSXX_R2_001.fastq.gz exome_fastq_merged/ECV2-31-plasma171215_tumor_tagseq-medexome_R2_001.fastq.gz &
# cp ECV2-31-plasma180202_tumor_takara-tagseq_BHGTHHDSXX_R2_001.fastq.gz exome_fastq_merged/ECV2-31-plasma180202_tumor_tagseq-medexome_R2_001.fastq.gz &
# cp ECV2-35-biopsi-G1_tumor_takara-tagseq_BHGTHHDSXX_R2_001.fastq.gz    exome_fastq_merged/ECV2-35-biopsi-G1_tumor_tagseq-medexome_R2_001.fastq.gz &
# cp ECV2-35-biopsi-G2_tumor_takara-tagseq_BHGTHHDSXX_R2_001.fastq.gz    exome_fastq_merged/ECV2-35-biopsi-G2_tumor_tagseq-medexome_R2_001.fastq.gz &
# cp ECV2-35-plasma180102_tumor_takara-tagseq_BHGTHHDSXX_R2_001.fastq.gz exome_fastq_merged/ECV2-35-plasma180102_tumor_tagseq-medexome_R2_001.fastq.gz &
# cp ECV2-35-plasma180316_tumor_takara-tagseq_BHGTHHDSXX_R2_001.fastq.gz exome_fastq_merged/ECV2-35-plasma180316_tumor_tagseq-medexome_R2_001.fastq.gz &
# cp ECV2-35-plasma180601_tumor_takara-tagseq_BHGTHHDSXX_R2_001.fastq.gz exome_fastq_merged/ECV2-35-plasma180601_tumor_tagseq-medexome_R2_001.fastq.gz &
# cp ECV2-35-recidiv_takara-tagseq_BHGTHHDSXX_R2_001.fastq.gz            exome_fastq_merged/ECV2-35-recidiv_tagseq-medexome_R2_001.fastq.gz &
# cp ECV2-4-biopsi-H1_tumor_takara-tagseq_BHGTHHDSXX_R2_001.fastq.gz     exome_fastq_merged/ECV2-4-biopsi-H1_tumor_tagseq-medexome_R2_001.fastq.gz &
# cp ECV2-4-biopsi-H2_tumor_takara-tagseq_BHGTHHDSXX_R2_001.fastq.gz     exome_fastq_merged/ECV2-4-biopsi-H2_tumor_tagseq-medexome_R2_001.fastq.gz &
# cp ECV2-4-blod_normal_takara-tagseq_BHGTHHDSXX_R2_001.fastq.gz         exome_fastq_merged/ECV2-4-blod_normal_tagseq-medexome_R2_001.fastq.gz &
# cp ECV2-4-op-1401_tumor_takara-tagseq_BHGTHHDSXX_R2_001.fastq.gz       exome_fastq_merged/ECV2-4-op-1401_tumor_tagseq-medexome_R2_001.fastq.gz &
# cp ECV2-4-op-1801_tumor_takara-tagseq_BHGTHHDSXX_R2_001.fastq.gz       exome_fastq_merged/ECV2-4-op-1801_tumor_tagseq-medexome_R2_001.fastq.gz &
# cp ECV2-4-plasma170712_tumor_takara-tagseq_BHGTHHDSXX_R2_001.fastq.gz  exome_fastq_merged/ECV2-4-plasma170712_tumor_tagseq-medexome_R2_001.fastq.gz &
# cp ECV2-4-plasma180205_tumor_takara-tagseq_BHGTHHDSXX_R2_001.fastq.gz  exome_fastq_merged/ECV2-4-plasma180205_tumor_tagseq-medexome_R2_001.fastq.gz &
# cp ECV2-8-biopsi-I1_tumor_takara-tagseq_BHGTHHDSXX_R2_001.fastq.gz     exome_fastq_merged/ECV2-8-biopsi-I1_tumor_tagseq-medexome_R2_001.fastq.gz &
# cp ECV2-8-biopsi-I2_tumor_takara-tagseq_BHGTHHDSXX_R2_001.fastq.gz     exome_fastq_merged/ECV2-8-biopsi-I2_tumor_tagseq-medexome_R2_001.fastq.gz &
# cp ECV2-8-blod_normal_takara-tagseq_BHGTHHDSXX_R2_001.fastq.gz         exome_fastq_merged/ECV2-8-blod_normal_tagseq-medexome_R2_001.fastq.gz &
# cp ECV2-8-op-01_tumor_takara-tagseq_BHGTHHDSXX_R2_001.fastq.gz         exome_fastq_merged/ECV2-8-op-01_tumor_tagseq-medexome_R2_001.fastq.gz &
# cp ECV2-8-op-02_tumor_takara-tagseq_BHGTHHDSXX_R2_001.fastq.gz         exome_fastq_merged/ECV2-8-op-02_tumor_tagseq-medexome_R2_001.fastq.gz &
# cp ECV2-8-plasma170522_tumor_takara-tagseq_BHGTHHDSXX_R2_001.fastq.gz  exome_fastq_merged/ECV2-8-plasma170522_tumor_tagseq-medexome_R2_001.fastq.gz &
# cp ECV2-8-plasma170728_tumor_takara-tagseq_BHGTHHDSXX_R2_001.fastq.gz  exome_fastq_merged/ECV2-8-plasma170728_tumor_tagseq-medexome_R2_001.fastq.gz &
# cp ECV2-8-plasma180806_tumor_takara-tagseq_BHGTHHDSXX_R2_001.fastq.gz  exome_fastq_merged/ECV2-8-plasma180806_tumor_tagseq-medexome_R2_001.fastq.gz &
# cp PC1-10-blod_normal_takara-tagseq_BHGTHHDSXX_R2_001.fastq.gz         exome_fastq_merged/PC1-10-blod_normal_tagseq-medexome_R2_001.fastq.gz &
# cp PC1-10-op-P1_tumor_takara-tagseq_BHGTHHDSXX_R2_001.fastq.gz         exome_fastq_merged/PC1-10-op-P1_tumor_tagseq-medexome_R2_001.fastq.gz &
# cp PC1-10-op-P2_tumor_takara-tagseq_BHGTHHDSXX_R2_001.fastq.gz         exome_fastq_merged/PC1-10-op-P2_tumor_tagseq-medexome_R2_001.fastq.gz &
# cp PC1-10-op-P3_tumor_takara-tagseq_BHGTHHDSXX_R2_001.fastq.gz         exome_fastq_merged/PC1-10-op-P3_tumor_tagseq-medexome_R2_001.fastq.gz &
# cp PC1-10-op-P4_tumor_takara-tagseq_BHGTHHDSXX_R2_001.fastq.gz         exome_fastq_merged/PC1-10-op-P4_tumor_tagseq-medexome_R2_001.fastq.gz &
# cp PC1-10-plasma170928_tumor_takara-tagseq_BHGTHHDSXX_R2_001.fastq.gz  exome_fastq_merged/PC1-10-plasma170928_tumor_tagseq-medexome_R2_001.fastq.gz &
# cp PC1-14-blod_normal_takara-tagseq_BHGTHHDSXX_R2_001.fastq.gz         exome_fastq_merged/PC1-14-blod_normal_tagseq-medexome_R2_001.fastq.gz &
# cp PC1-14-op-P5_tumor_takara-tagseq_BHGTHHDSXX_R2_001.fastq.gz         exome_fastq_merged/PC1-14-op-P5_tumor_tagseq-medexome_R2_001.fastq.gz &
# cp PC1-14-op-P6_tumor_takara-tagseq_BHGTHHDSXX_R2_001.fastq.gz         exome_fastq_merged/PC1-14-op-P6_tumor_tagseq-medexome_R2_001.fastq.gz &
# cp PC1-14-op-P7_tumor_takara-tagseq_BHGTHHDSXX_R2_001.fastq.gz         exome_fastq_merged/PC1-14-op-P7_tumor_tagseq-medexome_R2_001.fastq.gz &
# cp PC1-14-op-P8_tumor_takara-tagseq_BHGTHHDSXX_R2_001.fastq.gz         exome_fastq_merged/PC1-14-op-P8_tumor_tagseq-medexome_R2_001.fastq.gz &
# cp PC1-14-plasma180111_tumor_takara-tagseq_BHGTHHDSXX_R2_001.fastq.gz  exome_fastq_merged/PC1-14-plasma180111_tumor_tagseq-medexome_R2_001.fastq.gz &
# cp PC1-14-plasma180907_tumor_takara-tagseq_BHGTHHDSXX_R2_001.fastq.gz  exome_fastq_merged/PC1-14-plasma180907_tumor_tagseq-medexome_R2_001.fastq.gz &
# cp PC1-14-recidiv_takara-tagseq_BHGTHHDSXX_R2_001.fastq.gz             exome_fastq_merged/PC1-14-recidiv_tagseq-medexome_R2_001.fastq.gz &
# cp PC1-18-blod_normal_takara-tagseq_BHGTHHDSXX_R2_001.fastq.gz         exome_fastq_merged/PC1-18-blod_normal_tagseq-medexome_R2_001.fastq.gz &
# cp PC1-18-op-P10_tumor_takara-tagseq_BHGTHHDSXX_R2_001.fastq.gz        exome_fastq_merged/PC1-18-op-P10_tumor_tagseq-medexome_R2_001.fastq.gz &
# cp PC1-18-op-P11_tumor_takara-tagseq_BHGTHHDSXX_R2_001.fastq.gz        exome_fastq_merged/PC1-18-op-P11_tumor_tagseq-medexome_R2_001.fastq.gz &
# cp PC1-18-op-P9_tumor_takara-tagseq_BHGTHHDSXX_R2_001.fastq.gz         exome_fastq_merged/PC1-18-op-P9_tumor_tagseq-medexome_R2_001.fastq.gz &
# cp PC1-18-plasma180219_tumor_takara-tagseq_BHGTHHDSXX_R2_001.fastq.gz  exome_fastq_merged/PC1-18-plasma180219_tumor_tagseq-medexome_R2_001.fastq.gz &

# Sbatch
zcat ECV2-29-biopsi-C1_tumor_tagseq-medexome_HJ7KLBGX91_R1_001.fastq.gz  ECV2-29-biopsi-C1_tumor_takara-tagseq_BHGTHHDSXX_R1_001.fastq.gz 	| gzip -c > exome_fastq_merged/ECV2-29-biopsi-C1_tumor_tagseq-medexome_R1_001.fastq.gz &
zcat ECV2-29-biopsi-C2_tumor_tagseq-medexome_HJ7KLBGX92_R1_001.fastq.gz  ECV2-29-biopsi-C2_tumor_takara-tagseq_BHGTHHDSXX_R1_001.fastq.gz	| gzip -c > exome_fastq_merged/ECV2-29-biopsi-C2_tumor_tagseq-medexome_R1_001.fastq.gz &
zcat ECV2-29-op-A1_tumor_tagseq-medexome_HJ7KLBGX9_R1_001.fastq.gz       ECV2-29-op-A1_tumor_takara-tagseq_BHGTHHDSXX_R1_001.fastq.gz		| gzip -c > exome_fastq_merged/ECV2-29-op-A1_tumor_tagseq-medexome_R1_001.fastq.gz &
zcat ECV2-29-op-A2_tumor_tagseq-medexome_HJ7KLBGX9_R1_001.fastq.gz       ECV2-29-op-A2_tumor_takara-tagseq_BHGTHHDSXX_R1_001.fastq.gz		| gzip -c > exome_fastq_merged/ECV2-29-op-A2_tumor_tagseq-medexome_R1_001.fastq.gz &
zcat ECV2-29-op-B1_tumor_tagseq-medexome_HJ7KLBGX9_R1_001.fastq.gz       ECV2-29-op-B1_tumor_takara-tagseq_BHGTHHDSXX_R1_001.fastq.gz		| gzip -c > exome_fastq_merged/ECV2-29-op-B1_tumor_tagseq-medexome_R1_001.fastq.gz &
zcat ECV2-29-op-B2_tumor_tagseq-medexome_HJ7KLBGX9_R1_001.fastq.gz       ECV2-29-op-B2_tumor_takara-tagseq_BHGTHHDSXX_R1_001.fastq.gz		| gzip -c > exome_fastq_merged/ECV2-29-op-B2_tumor_tagseq-medexome_R1_001.fastq.gz &
zcat ECV2-31-op-D2_tumor_tagseq-medexome_HC3NJBGX9_R1_001.fastq.gz       ECV2-31-op-D2_tumor_takara-tagseq_BHGTHHDSXX_R1_001.fastq.gz		| gzip -c > exome_fastq_merged/ECV2-31-op-D2_tumor_tagseq-medexome_R1_001.fastq.gz &
zcat ECV2-35-blod_normal_tagseq-medexome_HC3NJBGX9_R1_001.fastq.gz       ECV2-35-blod_normal_takara-tagseq_BHGTHHDSXX_R1_001.fastq.gz		| gzip -c > exome_fastq_merged/ECV2-35-blod_normal_tagseq-medexome_R1_001.fastq.gz &
zcat ECV2-35-op-01_tumor_tagseq-medexome_HJ7KLBGX9_R1_001.fastq.gz       ECV2-35-op-01_tumor_takara-tagseq_BHGTHHDSXX_R1_001.fastq.gz 		| gzip -c > exome_fastq_merged/ECV2-35-op-01_tumor_tagseq-medexome_R1_001.fastq.gz &
zcat ECV2-35-op-02_tumor_tagseq-medexome_HJ7KLBGX9_R1_001.fastq.gz       ECV2-35-op-02_tumor_takara-tagseq_BHGTHHDSXX_R1_001.fastq.gz		| gzip -c > exome_fastq_merged/ECV2-35-op-02_tumor_tagseq-medexome_R1_001.fastq.gz &
zcat ECV2-35-op-03_tumor_tagseq-medexome_HJ7KLBGX9_R1_001.fastq.gz       ECV2-35-op-03_tumor_takara-tagseq_BHGTHHDSXX_R1_001.fastq.gz		| gzip -c > exome_fastq_merged/ECV2-35-op-03_tumor_tagseq-medexome_R1_001.fastq.gz &
zcat ECV2-35-op-04_tumor_tagseq-medexome_HJ7KLBGX9_R1_001.fastq.gz       ECV2-35-op-04_tumor_takara-tagseq_BHGTHHDSXX_R1_001.fastq.gz		| gzip -c > exome_fastq_merged/ECV2-35-op-04_tumor_tagseq-medexome_R1_001.fastq.gz &
zcat ECV2-35-op-06_tumor_tagseq-medexome_HJ7KLBGX9_R1_001.fastq.gz       ECV2-35-op-06_tumor_takara-tagseq_BHGTHHDSXX_R1_001.fastq.gz		| gzip -c > exome_fastq_merged/ECV2-35-op-06_tumor_tagseq-medexome_R1_001.fastq.gz &
zcat ECV2-35-op-07_tumor_tagseq-medexome_HJ7KLBGX9_R1_001.fastq.gz       ECV2-35-op-07_tumor_takara-tagseq_BHGTHHDSXX_R1_001.fastq.gz		| gzip -c > exome_fastq_merged/ECV2-35-op-07_tumor_tagseq-medexome_R1_001.fastq.gz &
zcat ECV2-4-op-1802_tumor_tagseq-medexome_HC3NJBGX9_R1_001.fastq.gz      ECV2-4-op-1802_tumor_takara-tagseq_BHGTHHDSXX_R1_001.fastq.gz		| gzip -c > exome_fastq_merged/ECV2-4-op-1802_tumor_tagseq-medexome_R1_001.fastq.gz &
zcat ECV2-4-plasma170503_tumor_tagseq-medexome_HC3NJBGX9_R1_001.fastq.gz ECV2-4-plasma170503_tumor_takara-tagseq_BHGTHHDSXX_R1_001.fastq.gz	| gzip -c > exome_fastq_merged/ECV2-4-plasma170503_tumor_tagseq-medexome_R1_001.fastq.gz &
zcat ECV2-8-op-03_tumor_tagseq-medexome_HC3NJBGX9_R1_001.fastq.gz        ECV2-8-op-03_tumor_takara-tagseq_BHGTHHDSXX_R1_001.fastq.gz	    | gzip -c > exome_fastq_merged/ECV2-8-op-03_tumor_tagseq-medexome_R1_001.fastq.gz &
zcat PC1-18-op-P12_tumor_tagseq-medexome_HC3NJBGX9_R1_001.fastq.gz       PC1-18-op-P12_tumor_takara-tagseq_BHGTHHDSXX_R1_001.fastq.gz		| gzip -c > exome_fastq_merged/PC1-18-op-P12_tumor_tagseq-medexome_R1_001.fastq.gz &

zcat ECV2-29-biopsi-C1_tumor_tagseq-medexome_HJ7KLBGX91_R2_001.fastq.gz  ECV2-29-biopsi-C1_tumor_takara-tagseq_BHGTHHDSXX_R2_001.fastq.gz 	| gzip -c > exome_fastq_merged/ECV2-29-biopsi-C1_tumor_tagseq-medexome_R2_001.fastq.gz &
zcat ECV2-29-biopsi-C2_tumor_tagseq-medexome_HJ7KLBGX92_R2_001.fastq.gz  ECV2-29-biopsi-C2_tumor_takara-tagseq_BHGTHHDSXX_R2_001.fastq.gz	| gzip -c > exome_fastq_merged/ECV2-29-biopsi-C2_tumor_tagseq-medexome_R2_001.fastq.gz &
zcat ECV2-29-op-A1_tumor_tagseq-medexome_HJ7KLBGX9_R2_001.fastq.gz       ECV2-29-op-A1_tumor_takara-tagseq_BHGTHHDSXX_R2_001.fastq.gz		| gzip -c > exome_fastq_merged/ECV2-29-op-A1_tumor_tagseq-medexome_R2_001.fastq.gz &
zcat ECV2-29-op-A2_tumor_tagseq-medexome_HJ7KLBGX9_R2_001.fastq.gz       ECV2-29-op-A2_tumor_takara-tagseq_BHGTHHDSXX_R2_001.fastq.gz		| gzip -c > exome_fastq_merged/ECV2-29-op-A2_tumor_tagseq-medexome_R2_001.fastq.gz &
zcat ECV2-29-op-B1_tumor_tagseq-medexome_HJ7KLBGX9_R2_001.fastq.gz       ECV2-29-op-B1_tumor_takara-tagseq_BHGTHHDSXX_R2_001.fastq.gz		| gzip -c > exome_fastq_merged/ECV2-29-op-B1_tumor_tagseq-medexome_R2_001.fastq.gz &
zcat ECV2-29-op-B2_tumor_tagseq-medexome_HJ7KLBGX9_R2_001.fastq.gz       ECV2-29-op-B2_tumor_takara-tagseq_BHGTHHDSXX_R2_001.fastq.gz		| gzip -c > exome_fastq_merged/ECV2-29-op-B2_tumor_tagseq-medexome_R2_001.fastq.gz &
sleep 18000

# On log in node
# zcat ECV2-31-op-D2_tumor_tagseq-medexome_HC3NJBGX9_R2_001.fastq.gz       ECV2-31-op-D2_tumor_takara-tagseq_BHGTHHDSXX_R2_001.fastq.gz		| gzip -c > exome_fastq_merged/ECV2-31-op-D2_tumor_tagseq-medexome_R2_001.fastq.gz &
# zcat ECV2-35-blod_normal_tagseq-medexome_HC3NJBGX9_R2_001.fastq.gz       ECV2-35-blod_normal_takara-tagseq_BHGTHHDSXX_R2_001.fastq.gz		| gzip -c > exome_fastq_merged/ECV2-35-blod_normal_tagseq-medexome_R2_001.fastq.gz &
# zcat ECV2-35-op-01_tumor_tagseq-medexome_HJ7KLBGX9_R2_001.fastq.gz       ECV2-35-op-01_tumor_takara-tagseq_BHGTHHDSXX_R2_001.fastq.gz 		| gzip -c > exome_fastq_merged/ECV2-35-op-01_tumor_tagseq-medexome_R2_001.fastq.gz &
# zcat ECV2-35-op-02_tumor_tagseq-medexome_HJ7KLBGX9_R2_001.fastq.gz       ECV2-35-op-02_tumor_takara-tagseq_BHGTHHDSXX_R2_001.fastq.gz		| gzip -c > exome_fastq_merged/ECV2-35-op-02_tumor_tagseq-medexome_R2_001.fastq.gz &
# zcat ECV2-35-op-03_tumor_tagseq-medexome_HJ7KLBGX9_R2_001.fastq.gz       ECV2-35-op-03_tumor_takara-tagseq_BHGTHHDSXX_R2_001.fastq.gz		| gzip -c > exome_fastq_merged/ECV2-35-op-03_tumor_tagseq-medexome_R2_001.fastq.gz &
# zcat ECV2-35-op-04_tumor_tagseq-medexome_HJ7KLBGX9_R2_001.fastq.gz       ECV2-35-op-04_tumor_takara-tagseq_BHGTHHDSXX_R2_001.fastq.gz		| gzip -c > exome_fastq_merged/ECV2-35-op-04_tumor_tagseq-medexome_R2_001.fastq.gz &
# zcat ECV2-35-op-06_tumor_tagseq-medexome_HJ7KLBGX9_R2_001.fastq.gz       ECV2-35-op-06_tumor_takara-tagseq_BHGTHHDSXX_R2_001.fastq.gz		| gzip -c > exome_fastq_merged/ECV2-35-op-06_tumor_tagseq-medexome_R2_001.fastq.gz &
# zcat ECV2-35-op-07_tumor_tagseq-medexome_HJ7KLBGX9_R2_001.fastq.gz       ECV2-35-op-07_tumor_takara-tagseq_BHGTHHDSXX_R2_001.fastq.gz		| gzip -c > exome_fastq_merged/ECV2-35-op-07_tumor_tagseq-medexome_R2_001.fastq.gz &
# zcat ECV2-4-op-1802_tumor_tagseq-medexome_HC3NJBGX9_R2_001.fastq.gz      ECV2-4-op-1802_tumor_takara-tagseq_BHGTHHDSXX_R2_001.fastq.gz		| gzip -c > exome_fastq_merged/ECV2-4-op-1802_tumor_tagseq-medexome_R2_001.fastq.gz &
# zcat ECV2-4-plasma170503_tumor_tagseq-medexome_HC3NJBGX9_R2_001.fastq.gz ECV2-4-plasma170503_tumor_takara-tagseq_BHGTHHDSXX_R2_001.fastq.gz	| gzip -c > exome_fastq_merged/ECV2-4-plasma170503_tumor_tagseq-medexome_R2_001.fastq.gz &
# zcat ECV2-8-op-03_tumor_tagseq-medexome_HC3NJBGX9_R2_001.fastq.gz        ECV2-8-op-03_tumor_takara-tagseq_BHGTHHDSXX_R2_001.fastq.gz	    | gzip -c > exome_fastq_merged/ECV2-8-op-03_tumor_tagseq-medexome_R2_001.fastq.gz &
# zcat PC1-18-op-P12_tumor_tagseq-medexome_HC3NJBGX9_R2_001.fastq.gz       PC1-18-op-P12_tumor_takara-tagseq_BHGTHHDSXX_R2_001.fastq.gz		| gzip -c > exome_fastq_merged/PC1-18-op-P12_tumor_tagseq-medexome_R2_001.fastq.gz &
