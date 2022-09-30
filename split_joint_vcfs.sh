### MATCHED TUMOR-NORMAL VCFS FOR VAF-LOGR PLOTS ###
# OUT="single_sample_vcfs"
# OUT="single_sample_vcfs_NormADzero"
# OUT="single_sample_vcfs_pon-filtered_strict"
# EXTENSION="_somatic_mutect2.vcf.gz"
# EXTENSION="_cns_varscan2_somatic.vcf.gz"
# EXTENSION="_cns_varscan2_somatic_NormADzero.vcf.gz"
# EXTENSION="_cns_varscan2_somatic_pon-filtered_strict.vcf.gz"

### PAIRED VCFS FOR VAF-VAF PLOTS ###
# OUT="paired_tumor_samples_vcfs"
# OUT_EXTENSION="_somatic_mutect2.vcf.gz"
# IN_EXTENSION="_somatic_mutect2_filtered.vcf"
# OUT="paired_tumor_samples_vcfs_NormADzero"
# OUT_EXTENSION="_somatic_varscan2_NormAD0.vcf.gz"
# IN_EXTENSION="_cns_varscan2_somatic_NormADzero.vcf.gz"
# OUT="paired_tumor_samples_vcfs_NormAD1"
# OUT_EXTENSION="_somatic_varscan2_NormAD1.vcf.gz"
# IN_EXTENSION="_cns_varscan2_somatic.vcf.gz"
OUT="paired_tumor_samples_vcfs_pon-filtered_strict"
OUT_EXTENSION="_somatic_varscan2_pon-filtered_strict.vcf.gz"
IN_EXTENSION="_cns_varscan2_somatic_pon-filtered_strict.vcf.gz"


###############################################
##### Split to tumor + matched normal vcf #####
###############################################

### Split joint vcf with multiple sample to single vcfs with data from one tumor sample and the matched normal sample
### Will produce and when trying to create a single-vcf for the blood sample and create an empty file
### Remember to remove empty vcf for blood samples subsequently

# for sample in `bcftools query -l ECV2-29$EXTENSION`; do
#   bcftools view -c1 -Oz -s $sample,ECV2-29-blod_normal_tagseq-medexome -o $OUT/$sample.vcf.gz ECV2-29$EXTENSION
# done
#
# for sample in `bcftools query -l ECV2-31$EXTENSION` ; do
#   bcftools view -c1 -Oz -s $sample,ECV2-31-blod_normal_tagseq-medexome -o $OUT/$sample.vcf.gz ECV2-31$EXTENSION
# done
#
# for sample in `bcftools query -l ECV2-35$EXTENSION` ; do
#   bcftools view -c1 -Oz -s $sample,ECV2-35-blod_normal_tagseq-medexome -o $OUT/$sample.vcf.gz ECV2-35$EXTENSION
# done
#
# for sample in `bcftools query -l ECV2-4$EXTENSION` ; do
#   bcftools view -c1 -Oz -s $sample,ECV2-4-blod_normal_tagseq-medexome -o $OUT/$sample.vcf.gz ECV2-4$EXTENSION
# done
#
# for sample in `bcftools query -l ECV2-8$EXTENSION` ; do
#   bcftools view -c1 -Oz -s $sample,ECV2-8-blod_normal_tagseq-medexome -o $OUT/$sample.vcf.gz ECV2-8$EXTENSION
# done
#
# for sample in `bcftools query -l PC1-10$EXTENSION` ; do
#   bcftools view -c1 -Oz -s $sample,PC1-10-blod_normal_tagseq-medexome -o $OUT/$sample.vcf.gz PC1-10$EXTENSION
# done
#
# for sample in `bcftools query -l PC1-14$EXTENSION` ; do
#   bcftools view -c1 -Oz -s $sample,PC1-14-blod_normal_tagseq-medexome -o $OUT/$sample.vcf.gz PC1-14$EXTENSION
# done
#
# for sample in `bcftools query -l PC1-18$EXTENSION` ; do
#   bcftools view -c1 -Oz -s $sample,PC1-18-blod_normal_tagseq-medexome -o $OUT/$sample.vcf.gz PC1-18$EXTENSION
# done

#######################################################################
##### Split to vcf with two paired tumor samples + matched normal #####
#######################################################################

##### BIOPSI vs BIOPSI #####

# bcftools view -c1 -Oz -s ECV2-29-biopsi-C1_tumor_tagseq-medexome-deep-seq,ECV2-29-biopsi-C2_tumor_tagseq-medexome-deep-seq,ECV2-29-blod_normal_tagseq-medexome \
# -o $OUT/ECV2-29-biopsi-C1_vs_ECV2-29-biopsi-C2$OUT_EXTENSION ECV2-29$IN_EXTENSION
#
# bcftools view -c1 -Oz -s ECV2-31-biopsi-F1_tumor_tagseq-medexome-deep-seq,ECV2-31-biopsi-F2_tumor_tagseq-medexome-deep-seq,ECV2-31-blod_normal_tagseq-medexome \
# -o $OUT/ECV2-31-biopsi-F1_vs_ECV2-31-biopsi-F2$OUT_EXTENSION ECV2-31$IN_EXTENSION
#
# bcftools view -c1 -Oz -s ECV2-35-biopsi-G1_tumor_tagseq-medexome-deep-seq,ECV2-35-biopsi-G2_tumor_tagseq-medexome-deep-seq,ECV2-35-blod_normal_tagseq-medexome \
# -o $OUT/ECV2-35-biopsi-G1_vs_ECV2-35-biopsi-G2$OUT_EXTENSION ECV2-35$IN_EXTENSION
#
# bcftools view -c1 -Oz -s ECV2-4-biopsi-H1_tumor_tagseq-medexome-deep-seq,ECV2-4-biopsi-H2_tumor_tagseq-medexome-deep-seq,ECV2-4-blod_normal_tagseq-medexome \
# -o $OUT/ECV2-4-biopsi-H1_vs_ECV2-4-biopsi-H2$OUT_EXTENSION ECV2-4$IN_EXTENSION
#
# bcftools view -c1 -Oz -s ECV2-8-biopsi-I1_tumor_tagseq-medexome-deep-seq,ECV2-8-biopsi-I2_tumor_tagseq-medexome-deep-seq,ECV2-8-blod_normal_tagseq-medexome \
# -o $OUT/ECV2-8-biopsi-I1_vs_ECV2-8-biopsi-I2$OUT_EXTENSION ECV2-8$IN_EXTENSION
#
##### OPERATION vs OPERATION #####

# bcftools view -c1 -Oz -s ECV2-29-op-A1_tumor_tagseq-medexome-deep-seq,ECV2-29-op-A2_tumor_tagseq-medexome-deep-seq,ECV2-29-blod_normal_tagseq-medexome \
# -o $OUT/ECV2-29-op-A1_vs_ECV2-29-op-A2$OUT_EXTENSION ECV2-29$IN_EXTENSION
# bcftools view -c1 -Oz -s ECV2-29-op-A1_tumor_tagseq-medexome-deep-seq,ECV2-29-op-B1_tumor_tagseq-medexome-deep-seq,ECV2-29-blod_normal_tagseq-medexome \
# -o $OUT/ECV2-29-op-A1_vs_ECV2-29-op-B1$OUT_EXTENSION ECV2-29$IN_EXTENSION
# bcftools view -c1 -Oz -s ECV2-29-op-A1_tumor_tagseq-medexome-deep-seq,ECV2-29-op-B2_tumor_tagseq-medexome-deep-seq,ECV2-29-blod_normal_tagseq-medexome \
# -o $OUT/ECV2-29-op-A1_vs_ECV2-29-op-B2$OUT_EXTENSION ECV2-29$IN_EXTENSION
# bcftools view -c1 -Oz -s ECV2-29-op-A2_tumor_tagseq-medexome-deep-seq,ECV2-29-op-B1_tumor_tagseq-medexome-deep-seq,ECV2-29-blod_normal_tagseq-medexome \
# -o $OUT/ECV2-29-op-A2_vs_ECV2-29-op-B1$OUT_EXTENSION ECV2-29$IN_EXTENSION
# bcftools view -c1 -Oz -s ECV2-29-op-A2_tumor_tagseq-medexome-deep-seq,ECV2-29-op-B2_tumor_tagseq-medexome-deep-seq,ECV2-29-blod_normal_tagseq-medexome \
# -o $OUT/ECV2-29-op-A2_vs_ECV2-29-op-B2$OUT_EXTENSION ECV2-29$IN_EXTENSION
# bcftools view -c1 -Oz -s ECV2-29-op-B1_tumor_tagseq-medexome-deep-seq,ECV2-29-op-B2_tumor_tagseq-medexome-deep-seq,ECV2-29-blod_normal_tagseq-medexome \
# -o $OUT/ECV2-29-op-B1_vs_ECV2-29-op-B2$OUT_EXTENSION ECV2-29$IN_EXTENSION
#
# bcftools view -c1 -Oz -s ECV2-31-op-D1_tumor_tagseq-medexome-deep-seq,ECV2-31-op-D2_tumor_tagseq-medexome-deep-seq,ECV2-31-blod_normal_tagseq-medexome \
# -o $OUT/ECV2-31-op-D1_vs_ECV2-31-op-D2$OUT_EXTENSION ECV2-31$IN_EXTENSION
# bcftools view -c1 -Oz -s ECV2-31-op-D1_tumor_tagseq-medexome-deep-seq,ECV2-31-op-E1_tumor_tagseq-medexome-deep-seq,ECV2-31-blod_normal_tagseq-medexome \
# -o $OUT/ECV2-31-op-D1_vs_ECV2-31-op-E1$OUT_EXTENSION ECV2-31$IN_EXTENSION
# bcftools view -c1 -Oz -s ECV2-31-op-D1_tumor_tagseq-medexome-deep-seq,ECV2-31-op-E2_tumor_tagseq-medexome-deep-seq,ECV2-31-blod_normal_tagseq-medexome \
# -o $OUT/ECV2-31-op-D1_vs_ECV2-31-op-E2$OUT_EXTENSION ECV2-31$IN_EXTENSION
# bcftools view -c1 -Oz -s ECV2-31-op-D2_tumor_tagseq-medexome-deep-seq,ECV2-31-op-E1_tumor_tagseq-medexome-deep-seq,ECV2-31-blod_normal_tagseq-medexome \
# -o $OUT/ECV2-31-op-D2_vs_ECV2-31-op-E1$OUT_EXTENSION ECV2-31$IN_EXTENSION
# bcftools view -c1 -Oz -s ECV2-31-op-D2_tumor_tagseq-medexome-deep-seq,ECV2-31-op-E2_tumor_tagseq-medexome-deep-seq,ECV2-31-blod_normal_tagseq-medexome \
# -o $OUT/ECV2-31-op-D2_vs_ECV2-31-op-E2$OUT_EXTENSION ECV2-31$IN_EXTENSION
# bcftools view -c1 -Oz -s ECV2-31-op-E1_tumor_tagseq-medexome-deep-seq,ECV2-31-op-E2_tumor_tagseq-medexome-deep-seq,ECV2-31-blod_normal_tagseq-medexome \
# -o $OUT/ECV2-31-op-E1_vs_ECV2-31-op-E2$OUT_EXTENSION ECV2-31$IN_EXTENSION
#
# bcftools view -c1 -Oz -s ECV2-35-op-01_tumor_tagseq-medexome-deep-seq,ECV2-35-op-02_tumor_tagseq-medexome-deep-seq,ECV2-35-blod_normal_tagseq-medexome \
# -o $OUT/ECV2-35-op-01_vs_ECV2-35-op-02$OUT_EXTENSION ECV2-35$IN_EXTENSION
# bcftools view -c1 -Oz -s ECV2-35-op-01_tumor_tagseq-medexome-deep-seq,ECV2-35-op-03_tumor_tagseq-medexome-deep-seq,ECV2-35-blod_normal_tagseq-medexome \
# -o $OUT/ECV2-35-op-01_vs_ECV2-35-op-03$OUT_EXTENSION ECV2-35$IN_EXTENSION
# bcftools view -c1 -Oz -s ECV2-35-op-01_tumor_tagseq-medexome-deep-seq,ECV2-35-op-04_tumor_tagseq-medexome-deep-seq,ECV2-35-blod_normal_tagseq-medexome \
# -o $OUT/ECV2-35-op-01_vs_ECV2-35-op-04$OUT_EXTENSION ECV2-35$IN_EXTENSION
# bcftools view -c1 -Oz -s ECV2-35-op-01_tumor_tagseq-medexome-deep-seq,ECV2-35-op-06_tumor_tagseq-medexome-deep-seq,ECV2-35-blod_normal_tagseq-medexome \
# -o $OUT/ECV2-35-op-01_vs_ECV2-35-op-06$OUT_EXTENSION ECV2-35$IN_EXTENSION
# bcftools view -c1 -Oz -s ECV2-35-op-01_tumor_tagseq-medexome-deep-seq,ECV2-35-op-07_tumor_tagseq-medexome-deep-seq,ECV2-35-blod_normal_tagseq-medexome \
# -o $OUT/ECV2-35-op-01_vs_ECV2-35-op-07$OUT_EXTENSION ECV2-35$IN_EXTENSION
# bcftools view -c1 -Oz -s ECV2-35-op-02_tumor_tagseq-medexome-deep-seq,ECV2-35-op-03_tumor_tagseq-medexome-deep-seq,ECV2-35-blod_normal_tagseq-medexome \
# -o $OUT/ECV2-35-op-02_vs_ECV2-35-op-03$OUT_EXTENSION ECV2-35$IN_EXTENSION
# bcftools view -c1 -Oz -s ECV2-35-op-02_tumor_tagseq-medexome-deep-seq,ECV2-35-op-04_tumor_tagseq-medexome-deep-seq,ECV2-35-blod_normal_tagseq-medexome \
# -o $OUT/ECV2-35-op-02_vs_ECV2-35-op-04$OUT_EXTENSION ECV2-35$IN_EXTENSION
# bcftools view -c1 -Oz -s ECV2-35-op-02_tumor_tagseq-medexome-deep-seq,ECV2-35-op-06_tumor_tagseq-medexome-deep-seq,ECV2-35-blod_normal_tagseq-medexome \
# -o $OUT/ECV2-35-op-02_vs_ECV2-35-op-06$OUT_EXTENSION ECV2-35$IN_EXTENSION
# bcftools view -c1 -Oz -s ECV2-35-op-02_tumor_tagseq-medexome-deep-seq,ECV2-35-op-07_tumor_tagseq-medexome-deep-seq,ECV2-35-blod_normal_tagseq-medexome \
# -o $OUT/ECV2-35-op-02_vs_ECV2-35-op-07$OUT_EXTENSION ECV2-35$IN_EXTENSION
# bcftools view -c1 -Oz -s ECV2-35-op-03_tumor_tagseq-medexome-deep-seq,ECV2-35-op-04_tumor_tagseq-medexome-deep-seq,ECV2-35-blod_normal_tagseq-medexome \
# -o $OUT/ECV2-35-op-03_vs_ECV2-35-op-04$OUT_EXTENSION ECV2-35$IN_EXTENSION
# bcftools view -c1 -Oz -s ECV2-35-op-03_tumor_tagseq-medexome-deep-seq,ECV2-35-op-06_tumor_tagseq-medexome-deep-seq,ECV2-35-blod_normal_tagseq-medexome \
# -o $OUT/ECV2-35-op-03_vs_ECV2-35-op-06$OUT_EXTENSION ECV2-35$IN_EXTENSION
# bcftools view -c1 -Oz -s ECV2-35-op-03_tumor_tagseq-medexome-deep-seq,ECV2-35-op-07_tumor_tagseq-medexome-deep-seq,ECV2-35-blod_normal_tagseq-medexome \
# -o $OUT/ECV2-35-op-03_vs_ECV2-35-op-07$OUT_EXTENSION ECV2-35$IN_EXTENSION
# bcftools view -c1 -Oz -s ECV2-35-op-04_tumor_tagseq-medexome-deep-seq,ECV2-35-op-06_tumor_tagseq-medexome-deep-seq,ECV2-35-blod_normal_tagseq-medexome \
# -o $OUT/ECV2-35-op-04_vs_ECV2-35-op-06$OUT_EXTENSION ECV2-35$IN_EXTENSION
# bcftools view -c1 -Oz -s ECV2-35-op-04_tumor_tagseq-medexome-deep-seq,ECV2-35-op-07_tumor_tagseq-medexome-deep-seq,ECV2-35-blod_normal_tagseq-medexome \
# -o $OUT/ECV2-35-op-04_vs_ECV2-35-op-07$OUT_EXTENSION ECV2-35$IN_EXTENSION
# bcftools view -c1 -Oz -s ECV2-35-op-06_tumor_tagseq-medexome-deep-seq,ECV2-35-op-07_tumor_tagseq-medexome-deep-seq,ECV2-35-blod_normal_tagseq-medexome \
# -o $OUT/ECV2-35-op-06_vs_ECV2-35-op-07$OUT_EXTENSION ECV2-35$IN_EXTENSION
#
# bcftools view -c1 -Oz -s ECV2-4-op-1401_tumor_tagseq-medexome-deep-seq,ECV2-4-op-1801_tumor_tagseq-medexome-deep-seq,ECV2-4-blod_normal_tagseq-medexome \
# -o $OUT/ECV2-4-op-1401_vs_ECV2-4-op-1801$OUT_EXTENSION ECV2-4$IN_EXTENSION
# bcftools view -c1 -Oz -s ECV2-4-op-1401_tumor_tagseq-medexome-deep-seq,ECV2-4-op-1802_tumor_tagseq-medexome-deep-seq,ECV2-4-blod_normal_tagseq-medexome \
# -o $OUT/ECV2-4-op-1401_vs_ECV2-4-op-1802$OUT_EXTENSION ECV2-4$IN_EXTENSION
# bcftools view -c1 -Oz -s ECV2-4-op-1801_tumor_tagseq-medexome-deep-seq,ECV2-4-op-1802_tumor_tagseq-medexome-deep-seq,ECV2-4-blod_normal_tagseq-medexome \
# -o $OUT/ECV2-4-op-1801_vs_ECV2-4-op-1802$OUT_EXTENSION ECV2-4$IN_EXTENSION
#
# bcftools view -c1 -Oz -s ECV2-8-op-01_tumor_tagseq-medexome-deep-seq,ECV2-8-op-02_tumor_tagseq-medexome-deep-seq,ECV2-8-blod_normal_tagseq-medexome \
# -o $OUT/ECV2-8-op-01_vs_ECV2-8-op-02$OUT_EXTENSION ECV2-8$IN_EXTENSION
# bcftools view -c1 -Oz -s ECV2-8-op-01_tumor_tagseq-medexome-deep-seq,ECV2-8-op-03_tumor_tagseq-medexome-deep-seq,ECV2-8-blod_normal_tagseq-medexome \
# -o $OUT/ECV2-8-op-01_vs_ECV2-8-op-03$OUT_EXTENSION ECV2-8$IN_EXTENSION
# bcftools view -c1 -Oz -s ECV2-8-op-02_tumor_tagseq-medexome-deep-seq,ECV2-8-op-03_tumor_tagseq-medexome-deep-seq,ECV2-8-blod_normal_tagseq-medexome \
# -o $OUT/ECV2-8-op-02_vs_ECV2-8-op-03$OUT_EXTENSION ECV2-8$IN_EXTENSION
#
# bcftools view -c1 -Oz -s PC1-10-op-P1_tumor_tagseq-medexome-deep-seq,PC1-10-op-P2_tumor_tagseq-medexome-deep-seq,PC1-10-blod_normal_tagseq-medexome \
# -o $OUT/PC1-10-op-P1_vs_PC1-10-op-P2$OUT_EXTENSION PC1-10$IN_EXTENSION
# bcftools view -c1 -Oz -s PC1-10-op-P1_tumor_tagseq-medexome-deep-seq,PC1-10-op-P3_tumor_tagseq-medexome-deep-seq,PC1-10-blod_normal_tagseq-medexome \
# -o $OUT/PC1-10-op-P1_vs_PC1-10-op-P3$OUT_EXTENSION PC1-10$IN_EXTENSION
# bcftools view -c1 -Oz -s PC1-10-op-P1_tumor_tagseq-medexome-deep-seq,PC1-10-op-P4_tumor_tagseq-medexome-deep-seq,PC1-10-blod_normal_tagseq-medexome \
# -o $OUT/PC1-10-op-P1_vs_PC1-10-op-P4$OUT_EXTENSION PC1-10$IN_EXTENSION
# bcftools view -c1 -Oz -s PC1-10-op-P2_tumor_tagseq-medexome-deep-seq,PC1-10-op-P3_tumor_tagseq-medexome-deep-seq,PC1-10-blod_normal_tagseq-medexome \
# -o $OUT/PC1-10-op-P2_vs_PC1-10-op-P3$OUT_EXTENSION PC1-10$IN_EXTENSION
# bcftools view -c1 -Oz -s PC1-10-op-P2_tumor_tagseq-medexome-deep-seq,PC1-10-op-P4_tumor_tagseq-medexome-deep-seq,PC1-10-blod_normal_tagseq-medexome \
# -o $OUT/PC1-10-op-P2_vs_PC1-10-op-P4$OUT_EXTENSION PC1-10$IN_EXTENSION
# bcftools view -c1 -Oz -s PC1-10-op-P3_tumor_tagseq-medexome-deep-seq,PC1-10-op-P4_tumor_tagseq-medexome-deep-seq,PC1-10-blod_normal_tagseq-medexome \
# -o $OUT/PC1-10-op-P3_vs_PC1-10-op-P4$OUT_EXTENSION PC1-10$IN_EXTENSION
#
# bcftools view -c1 -Oz -s PC1-14-op-P5_tumor_tagseq-medexome-deep-seq,PC1-14-op-P6_tumor_tagseq-medexome-deep-seq,PC1-14-blod_normal_tagseq-medexome \
# -o $OUT/PC1-14-op-P5_vs_PC1-14-op-P6$OUT_EXTENSION PC1-14$IN_EXTENSION
# bcftools view -c1 -Oz -s PC1-14-op-P5_tumor_tagseq-medexome-deep-seq,PC1-14-op-P7_tumor_tagseq-medexome-deep-seq,PC1-14-blod_normal_tagseq-medexome \
# -o $OUT/PC1-14-op-P5_vs_PC1-14-op-P7$OUT_EXTENSION PC1-14$IN_EXTENSION
# bcftools view -c1 -Oz -s PC1-14-op-P5_tumor_tagseq-medexome-deep-seq,PC1-14-op-P8_tumor_tagseq-medexome-deep-seq,PC1-14-blod_normal_tagseq-medexome \
# -o $OUT/PC1-14-op-P5_vs_PC1-14-op-P8$OUT_EXTENSION PC1-14$IN_EXTENSION
# bcftools view -c1 -Oz -s PC1-14-op-P6_tumor_tagseq-medexome-deep-seq,PC1-14-op-P7_tumor_tagseq-medexome-deep-seq,PC1-14-blod_normal_tagseq-medexome \
# -o $OUT/PC1-14-op-P6_vs_PC1-14-op-P7$OUT_EXTENSION PC1-14$IN_EXTENSION
# bcftools view -c1 -Oz -s PC1-14-op-P6_tumor_tagseq-medexome-deep-seq,PC1-14-op-P8_tumor_tagseq-medexome-deep-seq,PC1-14-blod_normal_tagseq-medexome \
# -o $OUT/PC1-14-op-P6_vs_PC1-14-op-P8$OUT_EXTENSION PC1-14$IN_EXTENSION
# bcftools view -c1 -Oz -s PC1-14-op-P7_tumor_tagseq-medexome-deep-seq,PC1-14-op-P8_tumor_tagseq-medexome-deep-seq,PC1-14-blod_normal_tagseq-medexome \
# -o $OUT/PC1-14-op-P7_vs_PC1-14-op-P8$OUT_EXTENSION PC1-14$IN_EXTENSION
#
# bcftools view -c1 -Oz -s PC1-18-op-P10_tumor_tagseq-medexome-deep-seq,PC1-18-op-P11_tumor_tagseq-medexome-deep-seq,PC1-18-blod_normal_tagseq-medexome \
# -o $OUT/PC1-18-op-P10_vs_PC1-18-op-P11$OUT_EXTENSION PC1-18$IN_EXTENSION
# bcftools view -c1 -Oz -s PC1-18-op-P10_tumor_tagseq-medexome-deep-seq,PC1-18-op-P12_tumor_tagseq-medexome-deep-seq,PC1-18-blod_normal_tagseq-medexome \
# -o $OUT/PC1-18-op-P10_vs_PC1-18-op-P12$OUT_EXTENSION PC1-18$IN_EXTENSION
# bcftools view -c1 -Oz -s PC1-18-op-P10_tumor_tagseq-medexome-deep-seq,PC1-18-op-P9_tumor_tagseq-medexome-deep-seq,PC1-18-blod_normal_tagseq-medexome \
# -o $OUT/PC1-18-op-P10_vs_PC1-18-op-P9$OUT_EXTENSION PC1-18$IN_EXTENSION
# bcftools view -c1 -Oz -s PC1-18-op-P11_tumor_tagseq-medexome-deep-seq,PC1-18-op-P12_tumor_tagseq-medexome-deep-seq,PC1-18-blod_normal_tagseq-medexome \
# -o $OUT/PC1-18-op-P11_vs_PC1-18-op-P12$OUT_EXTENSION PC1-18$IN_EXTENSION
# bcftools view -c1 -Oz -s PC1-18-op-P11_tumor_tagseq-medexome-deep-seq,PC1-18-op-P9_tumor_tagseq-medexome-deep-seq,PC1-18-blod_normal_tagseq-medexome \
# -o $OUT/PC1-18-op-P11_vs_PC1-18-op-P9$OUT_EXTENSION PC1-18$IN_EXTENSION
# bcftools view -c1 -Oz -s PC1-18-op-P12_tumor_tagseq-medexome-deep-seq,PC1-18-op-P9_tumor_tagseq-medexome-deep-seq,PC1-18-blod_normal_tagseq-medexome \
# -o $OUT/PC1-18-op-P12_vs_PC1-18-op-P9$OUT_EXTENSION PC1-18$IN_EXTENSION
#
##### BIOPSI vs OPERATION #####

# bcftools view -c1 -Oz -s ECV2-29-biopsi-C1_tumor_tagseq-medexome-deep-seq,ECV2-29-op-A1_tumor_tagseq-medexome-deep-seq,ECV2-29-blod_normal_tagseq-medexome \
# -o $OUT/ECV2-29-biopsi-C1_vs_ECV2-29-op-A1$OUT_EXTENSION ECV2-29$IN_EXTENSION
# bcftools view -c1 -Oz -s ECV2-29-biopsi-C1_tumor_tagseq-medexome-deep-seq,ECV2-29-op-A2_tumor_tagseq-medexome-deep-seq,ECV2-29-blod_normal_tagseq-medexome \
# -o $OUT/ECV2-29-biopsi-C1_vs_ECV2-29-op-A2$OUT_EXTENSION ECV2-29$IN_EXTENSION
# bcftools view -c1 -Oz -s ECV2-29-biopsi-C1_tumor_tagseq-medexome-deep-seq,ECV2-29-op-B1_tumor_tagseq-medexome-deep-seq,ECV2-29-blod_normal_tagseq-medexome \
# -o $OUT/ECV2-29-biopsi-C1_vs_ECV2-29-op-B1$OUT_EXTENSION ECV2-29$IN_EXTENSION
# bcftools view -c1 -Oz -s ECV2-29-biopsi-C1_tumor_tagseq-medexome-deep-seq,ECV2-29-op-B2_tumor_tagseq-medexome-deep-seq,ECV2-29-blod_normal_tagseq-medexome \
# -o $OUT/ECV2-29-biopsi-C1_vs_ECV2-29-op-B2$OUT_EXTENSION ECV2-29$IN_EXTENSION
#
# bcftools view -c1 -Oz -s ECV2-29-biopsi-C2_tumor_tagseq-medexome-deep-seq,ECV2-29-op-A1_tumor_tagseq-medexome-deep-seq,ECV2-29-blod_normal_tagseq-medexome \
# -o $OUT/ECV2-29-biopsi-C2_vs_ECV2-29-op-A1$OUT_EXTENSION ECV2-29$IN_EXTENSION
# bcftools view -c1 -Oz -s ECV2-29-biopsi-C2_tumor_tagseq-medexome-deep-seq,ECV2-29-op-A2_tumor_tagseq-medexome-deep-seq,ECV2-29-blod_normal_tagseq-medexome \
# -o $OUT/ECV2-29-biopsi-C2_vs_ECV2-29-op-A2$OUT_EXTENSION ECV2-29$IN_EXTENSION
# bcftools view -c1 -Oz -s ECV2-29-biopsi-C2_tumor_tagseq-medexome-deep-seq,ECV2-29-op-B1_tumor_tagseq-medexome-deep-seq,ECV2-29-blod_normal_tagseq-medexome \
# -o $OUT/ECV2-29-biopsi-C2_vs_ECV2-29-op-B1$OUT_EXTENSION ECV2-29$IN_EXTENSION
# bcftools view -c1 -Oz -s ECV2-29-biopsi-C2_tumor_tagseq-medexome-deep-seq,ECV2-29-op-B2_tumor_tagseq-medexome-deep-seq,ECV2-29-blod_normal_tagseq-medexome \
# -o $OUT/ECV2-29-biopsi-C2_vs_ECV2-29-op-B2$OUT_EXTENSION ECV2-29$IN_EXTENSION
#
# bcftools view -c1 -Oz -s ECV2-31-biopsi-F1_tumor_tagseq-medexome-deep-seq,ECV2-31-op-D1_tumor_tagseq-medexome-deep-seq,ECV2-31-blod_normal_tagseq-medexome \
# -o $OUT/ECV2-31-biopsi-F1_vs_ECV2-31-op-D1$OUT_EXTENSION ECV2-31$IN_EXTENSION
# bcftools view -c1 -Oz -s ECV2-31-biopsi-F1_tumor_tagseq-medexome-deep-seq,ECV2-31-op-D2_tumor_tagseq-medexome-deep-seq,ECV2-31-blod_normal_tagseq-medexome \
# -o $OUT/ECV2-31-biopsi-F1_vs_ECV2-31-op-D2$OUT_EXTENSION ECV2-31$IN_EXTENSION
# bcftools view -c1 -Oz -s ECV2-31-biopsi-F1_tumor_tagseq-medexome-deep-seq,ECV2-31-op-E1_tumor_tagseq-medexome-deep-seq,ECV2-31-blod_normal_tagseq-medexome \
# -o $OUT/ECV2-31-biopsi-F1_vs_ECV2-31-op-E1$OUT_EXTENSION ECV2-31$IN_EXTENSION
# bcftools view -c1 -Oz -s ECV2-31-biopsi-F1_tumor_tagseq-medexome-deep-seq,ECV2-31-op-E2_tumor_tagseq-medexome-deep-seq,ECV2-31-blod_normal_tagseq-medexome \
# -o $OUT/ECV2-31-biopsi-F1_vs_ECV2-31-op-E2$OUT_EXTENSION ECV2-31$IN_EXTENSION
#
# bcftools view -c1 -Oz -s ECV2-31-biopsi-F2_tumor_tagseq-medexome-deep-seq,ECV2-31-op-D1_tumor_tagseq-medexome-deep-seq,ECV2-31-blod_normal_tagseq-medexome \
# -o $OUT/ECV2-31-biopsi-F2_vs_ECV2-31-op-D1$OUT_EXTENSION ECV2-31$IN_EXTENSION
# bcftools view -c1 -Oz -s ECV2-31-biopsi-F2_tumor_tagseq-medexome-deep-seq,ECV2-31-op-D2_tumor_tagseq-medexome-deep-seq,ECV2-31-blod_normal_tagseq-medexome \
# -o $OUT/ECV2-31-biopsi-F2_vs_ECV2-31-op-D2$OUT_EXTENSION ECV2-31$IN_EXTENSION
# bcftools view -c1 -Oz -s ECV2-31-biopsi-F2_tumor_tagseq-medexome-deep-seq,ECV2-31-op-E1_tumor_tagseq-medexome-deep-seq,ECV2-31-blod_normal_tagseq-medexome \
# -o $OUT/ECV2-31-biopsi-F2_vs_ECV2-31-op-E1$OUT_EXTENSION ECV2-31$IN_EXTENSION
# bcftools view -c1 -Oz -s ECV2-31-biopsi-F2_tumor_tagseq-medexome-deep-seq,ECV2-31-op-E2_tumor_tagseq-medexome-deep-seq,ECV2-31-blod_normal_tagseq-medexome \
# -o $OUT/ECV2-31-biopsi-F2_vs_ECV2-31-op-E2$OUT_EXTENSION ECV2-31$IN_EXTENSION
#
# bcftools view -c1 -Oz -s ECV2-35-biopsi-G1_tumor_tagseq-medexome-deep-seq,ECV2-35-op-01_tumor_tagseq-medexome-deep-seq,ECV2-35-blod_normal_tagseq-medexome \
# -o $OUT/ECV2-35-biopsi-G1_vs_ECV2-35-op-01$OUT_EXTENSION ECV2-35$IN_EXTENSION
# bcftools view -c1 -Oz -s ECV2-35-biopsi-G1_tumor_tagseq-medexome-deep-seq,ECV2-35-op-02_tumor_tagseq-medexome-deep-seq,ECV2-35-blod_normal_tagseq-medexome \
# -o $OUT/ECV2-35-biopsi-G1_vs_ECV2-35-op-02$OUT_EXTENSION ECV2-35$IN_EXTENSION
# bcftools view -c1 -Oz -s ECV2-35-biopsi-G1_tumor_tagseq-medexome-deep-seq,ECV2-35-op-03_tumor_tagseq-medexome-deep-seq,ECV2-35-blod_normal_tagseq-medexome \
# -o $OUT/ECV2-35-biopsi-G1_vs_ECV2-35-op-03$OUT_EXTENSION ECV2-35$IN_EXTENSION
# bcftools view -c1 -Oz -s ECV2-35-biopsi-G1_tumor_tagseq-medexome-deep-seq,ECV2-35-op-04_tumor_tagseq-medexome-deep-seq,ECV2-35-blod_normal_tagseq-medexome \
# -o $OUT/ECV2-35-biopsi-G1_vs_ECV2-35-op-04$OUT_EXTENSION ECV2-35$IN_EXTENSION
# bcftools view -c1 -Oz -s ECV2-35-biopsi-G1_tumor_tagseq-medexome-deep-seq,ECV2-35-op-06_tumor_tagseq-medexome-deep-seq,ECV2-35-blod_normal_tagseq-medexome \
# -o $OUT/ECV2-35-biopsi-G1_vs_ECV2-35-op-06$OUT_EXTENSION ECV2-35$IN_EXTENSION
# bcftools view -c1 -Oz -s ECV2-35-biopsi-G1_tumor_tagseq-medexome-deep-seq,ECV2-35-op-07_tumor_tagseq-medexome-deep-seq,ECV2-35-blod_normal_tagseq-medexome \
# -o $OUT/ECV2-35-biopsi-G1_vs_ECV2-35-op-07$OUT_EXTENSION ECV2-35$IN_EXTENSION
#
# bcftools view -c1 -Oz -s ECV2-35-biopsi-G2_tumor_tagseq-medexome-deep-seq,ECV2-35-op-01_tumor_tagseq-medexome-deep-seq,ECV2-35-blod_normal_tagseq-medexome \
# -o $OUT/ECV2-35-biopsi-G2_vs_ECV2-35-op-01$OUT_EXTENSION ECV2-35$IN_EXTENSION
# bcftools view -c1 -Oz -s ECV2-35-biopsi-G2_tumor_tagseq-medexome-deep-seq,ECV2-35-op-02_tumor_tagseq-medexome-deep-seq,ECV2-35-blod_normal_tagseq-medexome \
# -o $OUT/ECV2-35-biopsi-G2_vs_ECV2-35-op-02$OUT_EXTENSION ECV2-35$IN_EXTENSION
# bcftools view -c1 -Oz -s ECV2-35-biopsi-G2_tumor_tagseq-medexome-deep-seq,ECV2-35-op-03_tumor_tagseq-medexome-deep-seq,ECV2-35-blod_normal_tagseq-medexome \
# -o $OUT/ECV2-35-biopsi-G2_vs_ECV2-35-op-03$OUT_EXTENSION ECV2-35$IN_EXTENSION
# bcftools view -c1 -Oz -s ECV2-35-biopsi-G2_tumor_tagseq-medexome-deep-seq,ECV2-35-op-04_tumor_tagseq-medexome-deep-seq,ECV2-35-blod_normal_tagseq-medexome \
# -o $OUT/ECV2-35-biopsi-G2_vs_ECV2-35-op-04$OUT_EXTENSION ECV2-35$IN_EXTENSION
# bcftools view -c1 -Oz -s ECV2-35-biopsi-G2_tumor_tagseq-medexome-deep-seq,ECV2-35-op-06_tumor_tagseq-medexome-deep-seq,ECV2-35-blod_normal_tagseq-medexome \
# -o $OUT/ECV2-35-biopsi-G2_vs_ECV2-35-op-06$OUT_EXTENSION ECV2-35$IN_EXTENSION
# bcftools view -c1 -Oz -s ECV2-35-biopsi-G2_tumor_tagseq-medexome-deep-seq,ECV2-35-op-07_tumor_tagseq-medexome-deep-seq,ECV2-35-blod_normal_tagseq-medexome \
# -o $OUT/ECV2-35-biopsi-G2_vs_ECV2-35-op-07$OUT_EXTENSION ECV2-35$IN_EXTENSION
#
# bcftools view -c1 -Oz -s ECV2-4-biopsi-H1_tumor_tagseq-medexome-deep-seq,ECV2-4-op-1401_tumor_tagseq-medexome-deep-seq,ECV2-4-blod_normal_tagseq-medexome \
# -o $OUT/ECV2-4-biopsi-H1_vs_ECV2-4-op-1401$OUT_EXTENSION ECV2-4$IN_EXTENSION
# bcftools view -c1 -Oz -s ECV2-4-biopsi-H1_tumor_tagseq-medexome-deep-seq,ECV2-4-op-1801_tumor_tagseq-medexome-deep-seq,ECV2-4-blod_normal_tagseq-medexome \
# -o $OUT/ECV2-4-biopsi-H1_vs_ECV2-4-op-1801$OUT_EXTENSION ECV2-4$IN_EXTENSION
# bcftools view -c1 -Oz -s ECV2-4-biopsi-H1_tumor_tagseq-medexome-deep-seq,ECV2-4-op-1802_tumor_tagseq-medexome-deep-seq,ECV2-4-blod_normal_tagseq-medexome \
# -o $OUT/ECV2-4-biopsi-H1_vs_ECV2-4-op-1802$OUT_EXTENSION ECV2-4$IN_EXTENSION
#
# bcftools view -c1 -Oz -s ECV2-4-biopsi-H2_tumor_tagseq-medexome-deep-seq,ECV2-4-op-1401_tumor_tagseq-medexome-deep-seq,ECV2-4-blod_normal_tagseq-medexome \
# -o $OUT/ECV2-4-biopsi-H2_vs_ECV2-4-op-1401$OUT_EXTENSION ECV2-4$IN_EXTENSION
# bcftools view -c1 -Oz -s ECV2-4-biopsi-H2_tumor_tagseq-medexome-deep-seq,ECV2-4-op-1801_tumor_tagseq-medexome-deep-seq,ECV2-4-blod_normal_tagseq-medexome \
# -o $OUT/ECV2-4-biopsi-H2_vs_ECV2-4-op-1801$OUT_EXTENSION ECV2-4$IN_EXTENSION
# bcftools view -c1 -Oz -s ECV2-4-biopsi-H2_tumor_tagseq-medexome-deep-seq,ECV2-4-op-1802_tumor_tagseq-medexome-deep-seq,ECV2-4-blod_normal_tagseq-medexome \
# -o $OUT/ECV2-4-biopsi-H2_vs_ECV2-4-op-1802$OUT_EXTENSION ECV2-4$IN_EXTENSION
#
# bcftools view -c1 -Oz -s ECV2-8-biopsi-I1_tumor_tagseq-medexome-deep-seq,ECV2-8-op-01_tumor_tagseq-medexome-deep-seq,ECV2-8-blod_normal_tagseq-medexome \
# -o $OUT/ECV2-8-biopsi-I1_vs_ECV2-8-op-01$OUT_EXTENSION ECV2-8$IN_EXTENSION
# bcftools view -c1 -Oz -s ECV2-8-biopsi-I1_tumor_tagseq-medexome-deep-seq,ECV2-8-op-02_tumor_tagseq-medexome-deep-seq,ECV2-8-blod_normal_tagseq-medexome \
# -o $OUT/ECV2-8-biopsi-I1_vs_ECV2-8-op-02$OUT_EXTENSION ECV2-8$IN_EXTENSION
# bcftools view -c1 -Oz -s ECV2-8-biopsi-I1_tumor_tagseq-medexome-deep-seq,ECV2-8-op-03_tumor_tagseq-medexome-deep-seq,ECV2-8-blod_normal_tagseq-medexome \
# -o $OUT/ECV2-8-biopsi-I1_vs_ECV2-8-op-03$OUT_EXTENSION ECV2-8$IN_EXTENSION
#
# bcftools view -c1 -Oz -s ECV2-8-biopsi-I2_tumor_tagseq-medexome-deep-seq,ECV2-8-op-01_tumor_tagseq-medexome-deep-seq,ECV2-8-blod_normal_tagseq-medexome \
# -o $OUT/ECV2-8-biopsi-I2_vs_ECV2-8-op-01$OUT_EXTENSION ECV2-8$IN_EXTENSION
# bcftools view -c1 -Oz -s ECV2-8-biopsi-I2_tumor_tagseq-medexome-deep-seq,ECV2-8-op-02_tumor_tagseq-medexome-deep-seq,ECV2-8-blod_normal_tagseq-medexome \
# -o $OUT/ECV2-8-biopsi-I2_vs_ECV2-8-op-02$OUT_EXTENSION ECV2-8$IN_EXTENSION
# bcftools view -c1 -Oz -s ECV2-8-biopsi-I2_tumor_tagseq-medexome-deep-seq,ECV2-8-op-03_tumor_tagseq-medexome-deep-seq,ECV2-8-blod_normal_tagseq-medexome \
# -o $OUT/ECV2-8-biopsi-I2_vs_ECV2-8-op-03$OUT_EXTENSION ECV2-8$IN_EXTENSION
