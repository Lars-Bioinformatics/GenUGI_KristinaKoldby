### MATCHED TUMOR-NORMAL VCFS FOR VAF-LOGR PLOTS ###
# OUT="single_sample_vcfs"
# OUT="single_sample_vcfs_NormADzero"
OUT="single_sample_vcfs_pon-filtered_strict"
# EXTENSION="_somatic_mutect2.vcf.gz"
# EXTENSION="_cns_varscan2_somatic.vcf.gz"
# EXTENSION="_cns_varscan2_somatic_NormADzero.vcf.gz"
EXTENSION="_cns_varscan2_somatic_pon-filtered_strict.vcf.gz"

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
# OUT="paired_tumor_samples_vcfs_pon-filtered_strict"
# OUT_EXTENSION="_somatic_varscan2_pon-filtered_strict.vcf.gz"
# IN_EXTENSION="_cns_varscan2_somatic_pon-filtered_strict.vcf.gz"


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
for sample in `bcftools query -l ECV2-4$EXTENSION` ; do
  bcftools view -c1 -Oz -s $sample,ECV2-4-blod_normal_tagseq-medexome -o $OUT/$sample.vcf.gz ECV2-4$EXTENSION
done

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

##### OPERATION vs OPERATION #####

# bcftools view -c1 -Oz -s ECV2-29-opA-merged_tumor_tagseq-medexome-deep-seq,ECV2-29-opB-merged_tumor_tagseq-medexome-deep-seq,ECV2-29-blod_normal_tagseq-medexome \
# -o $OUT/ECV2-29-opA-merged_vs_ECV2-29-opB-merged$OUT_EXTENSION ECV2-29$IN_EXTENSION
#
# bcftools view -c1 -Oz -s ECV2-31-opD-merged_tumor_tagseq-medexome-deep-seq,ECV2-31-opE-merged_tumor_tagseq-medexome-deep-seq,ECV2-31-blod_normal_tagseq-medexome \
# -o $OUT/ECV2-31-opD-merged_vs_ECV2-31-opE-merged$OUT_EXTENSION ECV2-31$IN_EXTENSION
#
# bcftools view -c1 -Oz -s ECV2-35-op1to4-merged_tumor_tagseq-medexome-deep-seq,ECV2-35-op6to7-merged_tumor_tagseq-medexome-deep-seq,ECV2-35-blod_normal_tagseq-medexome \
# -o $OUT/ECV2-35-op1to4-merged_vs_ECV2-35-op6to7-merged$OUT_EXTENSION ECV2-35$IN_EXTENSION
#
# bcftools view -c1 -Oz -s ECV2-4-op18-merged_tumor_tagseq-medexome-deep-seq,ECV2-4-op-1401_tumor_tagseq-medexome-deep-seq,ECV2-4-blod_normal_tagseq-medexome \
# -o $OUT/ECV2-4-op18-merged_vs_ECV2-4-op-1401$OUT_EXTENSION ECV2-4$IN_EXTENSION


##### BIOPSI vs OPERATION #####

# bcftools view -c1 -Oz -s ECV2-29-biopsi-merged_tumor_tagseq-medexome-deep-seq,ECV2-29-opA-merged_tumor_tagseq-medexome-deep-seq,ECV2-29-blod_normal_tagseq-medexome \
# -o $OUT/ECV2-29-biopsi-merged_vs_ECV2-29-opA-merged$OUT_EXTENSION ECV2-29$IN_EXTENSION
#
# bcftools view -c1 -Oz -s ECV2-29-biopsi-merged_tumor_tagseq-medexome-deep-seq,ECV2-29-opB-merged_tumor_tagseq-medexome-deep-seq,ECV2-29-blod_normal_tagseq-medexome \
# -o $OUT/ECV2-29-biopsi-merged_vs_ECV2-29-opB-merged$OUT_EXTENSION ECV2-29$IN_EXTENSION
#
# bcftools view -c1 -Oz -s ECV2-31-biopsi-merged_tumor_tagseq-medexome-deep-seq,ECV2-31-opD-merged_tumor_tagseq-medexome-deep-seq,ECV2-31-blod_normal_tagseq-medexome \
# -o $OUT/ECV2-31-biopsi-merged_vs_ECV2-31-opD-merged$OUT_EXTENSION ECV2-31$IN_EXTENSION
#
# bcftools view -c1 -Oz -s ECV2-31-biopsi-merged_tumor_tagseq-medexome-deep-seq,ECV2-31-opE-merged_tumor_tagseq-medexome-deep-seq,ECV2-31-blod_normal_tagseq-medexome \
# -o $OUT/ECV2-31-biopsi-merged_vs_ECV2-31-opE-merged$OUT_EXTENSION ECV2-31$IN_EXTENSION
#
# bcftools view -c1 -Oz -s ECV2-35-biopsi-merged_tumor_tagseq-medexome-deep-seq,ECV2-35-op1to4-merged_tumor_tagseq-medexome-deep-seq,ECV2-35-blod_normal_tagseq-medexome \
# -o $OUT/ECV2-35-biopsi-merged_vs_ECV2-35-op1to4-merged$OUT_EXTENSION ECV2-35$IN_EXTENSION
#
# bcftools view -c1 -Oz -s ECV2-35-biopsi-merged_tumor_tagseq-medexome-deep-seq,ECV2-35-op6to7-merged_tumor_tagseq-medexome-deep-seq,ECV2-35-blod_normal_tagseq-medexome \
# -o $OUT/ECV2-35-biopsi-merged_vs_ECV2-35-op6to7-merged$OUT_EXTENSION ECV2-35$IN_EXTENSION
#
# bcftools view -c1 -Oz -s ECV2-4-biopsi-merged_tumor_tagseq-medexome-deep-seq,ECV2-4-op18-merged_tumor_tagseq-medexome-deep-seq,ECV2-4-blod_normal_tagseq-medexome \
# -o $OUT/ECV2-4-biopsi-merged_vs_ECV2-4-op18-merged$OUT_EXTENSION ECV2-4$IN_EXTENSION
#
# bcftools view -c1 -Oz -s ECV2-4-biopsi-merged_tumor_tagseq-medexome-deep-seq,ECV2-4-op-merged_tumor_tagseq-medexome-deep-seq,ECV2-4-blod_normal_tagseq-medexome \
# -o $OUT/ECV2-4-biopsi-merged_vs_ECV2-4-op-merged$OUT_EXTENSION ECV2-4$IN_EXTENSION
#
# bcftools view -c1 -Oz -s ECV2-8-biopsi-merged_tumor_tagseq-medexome-deep-seq,ECV2-8-op-merged_tumor_tagseq-medexome-deep-seq,ECV2-8-blod_normal_tagseq-medexome \
# -o $OUT/ECV2-8-biopsi-merged_vs_ECV2-8-op-merged$OUT_EXTENSION ECV2-8$IN_EXTENSION
