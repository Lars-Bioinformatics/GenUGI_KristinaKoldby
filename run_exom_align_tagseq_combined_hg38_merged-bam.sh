#conda activate snakemake

curDate=`date +"%F_%H-%M"`
log=log_exom_align_tagseq_combined_hg38_merged-bam_$curDate.log
echo 'Start:' $curDate > $log
# snakemake -s /work/sdukoldby/scripts/exom_align_tagseq_combined_hg38_merged-bam.smk -j 999 --cluster "sbatch --qos=long -A sdukoldby_fat --time 5-0" 2>&1 | tee -a $log
snakemake -s /work/sdukoldby/scripts/exom_align_tagseq_combined_hg38_merged-bam.smk -j 999 --cluster "sbatch -A sdukoldby_slim --time 24:00:00" 2>&1 | tee -a $log
echo 'Slut:' `date +%F_%H-%M` >> $log
