#conda activate snakemake

curDate=`date +"%F_%H-%M"`
log=log_wgs-align_qc_hg38_$curDate.log
echo 'Start:' $curDate > $log
snakemake -s /work/sdukoldby/scripts/wgs-align_qc_hg38.smk -j 999 --cluster "sbatch -A sdukoldby_slim --time 24:00:00" 2>&1 | tee -a $log
echo 'Slut:' `date +%F_%H-%M` >> $log
