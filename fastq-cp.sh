set -x

#cp ECV2-31-op-E1_tumor_takara-tagseq_BHGTHHDSXX_R1_001.fastq.gz exome_fastq_merged/ECV2-31-op-E1_tumor_tagseq-medexome_R1_001.fastq.gz

#cp ECV2-31-op-E1_tumor_takara-tagseq_BHGTHHDSXX_R2_001.fastq.gz exome_fastq_merged/ECV2-31-op-E1_tumor_tagseq-medexome_R2_001.fastq.gz

#cp ECV2-35-plasma180102_tumor_takara-tagseq_BHGTHHDSXX_R1_001.fastq.gz exome_fastq_merged/ECV2-35-plasma180102_tumor_tagseq-medexome_R1_001.fastq.gz

#cp ECV2-35-plasma180102_tumor_takara-tagseq_BHGTHHDSXX_R2_001.fastq.gz exome_fastq_merged/ECV2-35-plasma180102_tumor_tagseq-medexome_R2_001.fastq.gz

#cp ECV2-35-recidiv_takara-tagseq_BHGTHHDSXX_R1_001.fastq.gz exome_fastq_merged/ECV2-35-recidiv_tagseq-medexome_R1_001.fastq.gz

#cp ECV2-35-recidiv_takara-tagseq_BHGTHHDSXX_R2_001.fastq.gz exome_fastq_merged/ECV2-35-recidiv_tagseq-medexome_R2_001.fastq.gz

#cp ECV2-4-blod_normal_takara-tagseq_BHGTHHDSXX_R1_001.fastq.gz exome_fastq_merged/ECV2-4-blod_normal_tagseq-medexome_R1_001.fastq.gz

#cp ECV2-4-blod_normal_takara-tagseq_BHGTHHDSXX_R2_001.fastq.gz exome_fastq_merged/ECV2-4-blod_normal_tagseq-medexome_R2_001.fastq.gz

#cp ECV2-8-plasma170522_tumor_takara-tagseq_BHGTHHDSXX_R1_001.fastq.gz exome_fastq_merged/ECV2-8-plasma170522_tumor_tagseq-medexome_R1_001.fastq.gz

#cp ECV2-8-plasma170522_tumor_takara-tagseq_BHGTHHDSXX_R2_001.fastq.gz exome_fastq_merged/ECV2-8-plasma170522_tumor_tagseq-medexome_R2_001.fastq.gz

#cp PC1-14-recidiv_takara-tagseq_BHGTHHDSXX_R1_001.fastq.gz exome_fastq_merged/PC1-14-recidiv_tagseq-medexome_R1_001.fastq.gz

#cp PC1-14-recidiv_takara-tagseq_BHGTHHDSXX_R2_001.fastq.gz exome_fastq_merged/PC1-14-recidiv_tagseq-medexome_R2_001.fastq.gz

#cp PC1-18-op-P9_tumor_takara-tagseq_BHGTHHDSXX_R1_001.fastq.gz exome_fastq_merged/PC1-18-op-P9_tumor_tagseq-medexome_R1_001.fastq.gz

#cp PC1-18-op-P9_tumor_takara-tagseq_BHGTHHDSXX_R2_001.fastq.gz exome_fastq_merged/PC1-18-op-P9_tumor_tagseq-medexome_R2_001.fastq.gz

#zcat ECV2-29-op-A1_tumor_takara-tagseq_BHGTHHDSXX_R1_001.fastq.gz ECV2-29-op-A1_tumor_tagseq-medexome_HJ7KLBGX9_R1_001.fastq.gz | gzip -c > exome_fastq_merged/ECV2-29-op-A1_tumor_tagseq-medexome_R1_001.fastq.gz

#zcat ECV2-29-op-A1_tumor_takara-tagseq_BHGTHHDSXX_R2_001.fastq.gz ECV2-29-op-A1_tumor_tagseq-medexome_HJ7KLBGX9_R2_001.fastq.gz | gzip -c > exome_fastq_merged/ECV2-29-op-A1_tumor_tagseq-medexome_R2_001.fastq.gz

#zcat ECV2-29-op-A2_tumor_takara-tagseq_BHGTHHDSXX_R1_001.fastq.gz ECV2-29-op-A2_tumor_tagseq-medexome_HJ7KLBGX9_R1_001.fastq.gz | gzip -c > exome_fastq_merged/ECV2-29-op-A2_tumor_tagseq-medexome_R1_001.fastq.gz

#zcat ECV2-29-op-A2_tumor_takara-tagseq_BHGTHHDSXX_R2_001.fastq.gz ECV2-29-op-A2_tumor_tagseq-medexome_HJ7KLBGX9_R2_001.fastq.gz | gzip -c > exome_fastq_merged/ECV2-29-op-A2_tumor_tagseq-medexome_R2_001.fastq.gz

#zcat ECV2-29-op-B1_tumor_takara-tagseq_BHGTHHDSXX_R1_001.fastq.gz ECV2-29-op-B1_tumor_tagseq-medexome_HJ7KLBGX9_R1_001.fastq.gz | gzip -c > exome_fastq_merged/ECV2-29-op-B1_tumor_tagseq-medexome_R1_001.fastq.gz

#zcat ECV2-29-op-B1_tumor_takara-tagseq_BHGTHHDSXX_R2_001.fastq.gz ECV2-29-op-B1_tumor_tagseq-medexome_HJ7KLBGX9_R2_001.fastq.gz | gzip -c > exome_fastq_merged/ECV2-29-op-B1_tumor_tagseq-medexome_R2_001.fastq.gz

zcat ECV2-29-op-B2_tumor_takara-tagseq_BHGTHHDSXX_R1_001.fastq.gz ECV2-29-op-B2_tumor_tagseq-medexome_HJ7KLBGX9_R1_001.fastq.gz | gzip -c > exome_fastq_merged/ECV2-29-op-B2_tumor_tagseq-medexome_R1_001.fastq.gz &

zcat ECV2-29-op-B2_tumor_takara-tagseq_BHGTHHDSXX_R2_001.fastq.gz ECV2-29-op-B2_tumor_tagseq-medexome_HJ7KLBGX9_R2_001.fastq.gz | gzip -c > exome_fastq_merged/ECV2-29-op-B2_tumor_tagseq-medexome_R2_001.fastq.gz &

zcat ECV2-4-plasma170503_tumor_tagseq-medexome_HC3NJBGX9_R1_001.fastq.gz ECV2-4-plasma170503_tumor_takara-tagseq_BHGTHHDSXX_R1_001.fastq.gz | gzip -c > exome_fastq_merged/ECV2-4-plasma170503_tumor_tagseq-medexome_R1_001.fatsq.gz &

zcat ECV2-4-plasma170503_tumor_tagseq-medexome_HC3NJBGX9_R2_001.fastq.gz ECV2-4-plasma170503_tumor_takara-tagseq_BHGTHHDSXX_R2_001.fastq.gz | gzip -c > exome_fastq_merged/ECV2-4-plasma170503_tumor_tagseq-medexome_R2_001.fatsq.gz &