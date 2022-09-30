#conda activate snakemake

curDate=`date +"%F_%H-%M"`
log=log_exom_tagseq_connor_hg38_$curDate.log
echo 'Start:' $curDate > $log
snakemake -s /work/sdukoldby/scripts/exom_tagseq_connor_hg38.smk -j 999 --cluster "sbatch -A sdukoldby_fat --time 08:00:00" 2>&1 | tee -a $log
# snakemake -s /work/sdukoldby/scripts/exom_tagseq_connor_hg38.smk -j 999 --cluster "sbatch -A sdukoldby_fat --qos=long --time 2-0" 2>&1 | tee -a $log
echo 'Slut:' `date +%F_%H-%M` >> $log
