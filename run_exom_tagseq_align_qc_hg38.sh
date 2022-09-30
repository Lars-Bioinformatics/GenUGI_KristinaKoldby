#conda activate snakemake

curDate=`date +"%F_%H-%M"`
log=log_exom_tagseq_align_qc_hg38_$curDate.log
echo 'Start:' $curDate > $log
snakemake -s /work/sdukoldby/scripts/exom_tagseq_align_qc_hg38.smk 2>&1 | tee -a $log
echo 'Slut:' `date +%F_%H-%M` >> $log

#-j 999 --cluster "sbatch -A sdukoldby_slim --time 1-0"