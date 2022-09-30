curDate=`date +"%F_%H-%M"`
log=log_ngcgh_$curDate.log
echo "Start:" $curDate > $log
snakemake -s /work/sdukoldby/scripts/ngcgh.smk -j 999 --cluster "sbatch -A sdukoldby_slim --time 1-0" 2>&1 | tee -a $log
echo "Slut:" `date +%F_%H-%M` >> $log