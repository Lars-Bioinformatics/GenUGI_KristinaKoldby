cd /work/sduhumac/kristina/data/genova/plasma/


#LOGFILE="/scratch/sduhumac/kristina/data/genova/plasma/output/logfiles/G35-1_illumina-truseq-genome_merged_log_merge_bam_markduplicates.txt"
sbatch /gpfs/gss1/work/sduhumac/kristina/scripts/merge_bam_markduplicates.sh G35-1_illumina-truseq-genome_merged G35-1_illumina-truseq-genome_HTHCHBGXY G35-1_illumina-truseq-genome_HTJFVBGXY
#sbatch /gpfs/gss1/work/sduhumac/kristina/scripts/merge_bam_markduplicates.sh G35-2_illumina-truseq-genome_merged G35-2_illumina-truseq-genome_H3YVTBGX2 G35-2_illumina-truseq-genome_HTJW2BGXY
#sbatch /gpfs/gss1/work/sduhumac/kristina/scripts/merge_bam_markduplicates.sh G35-3_illumina-truseq-genome_HY35GBGXY
sbatch /gpfs/gss1/work/sduhumac/kristina/scripts/merge_bam_markduplicates.sh G35-4_illumina-truseq-genome_merged G35-4_illumina-truseq-genome_H3WWJBGX2 G35-4_illumina-truseq-genome_HTKG2BGXY
#sbatch /gpfs/gss1/work/sduhumac/kristina/scripts/merge_bam_markduplicates.sh G35-5_illumina-truseq-genome_HY22CBGXY 

#for i in *_R1.fastq.gz;
#do
#SAMPLEID=$(basename $i _R1.fastq.gz)
#echo $SAMPLEID
#if ! [[ $SAMPLEID == "G35-2"* ]]; then
#sbatch /work/sduhumac/kristina/scripts/merge_bam_markduplicates.sh $SAMPLEID
#echo $SAMPLEID
#fi
#done