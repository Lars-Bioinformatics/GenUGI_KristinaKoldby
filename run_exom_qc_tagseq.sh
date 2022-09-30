#conda activate snakemake

curDate=`date +"%F_%H-%M"`
log=log_exom_qc_tagseq_$curDate.log
echo 'Start:' $curDate > $log
snakemake -s /work/sdukoldby/scripts/exom_qc_tagseq.smk -j 999 --cluster "sbatch -A sdukoldby_slim --time 3:00:00" 2>&1 | tee -a $log
echo 'Slut:' `date +%F_%H-%M` >> $log

