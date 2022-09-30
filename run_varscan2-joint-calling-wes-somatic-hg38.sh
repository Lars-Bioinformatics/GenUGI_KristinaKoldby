curDate=`date +"%F_%H-%M"`
# log=log_varscan2-joint-calling-wes-somatic-hg38_$curDate.log
log=log_varscan2-joint-calling-deep-wes-somatic-hg38_$curDate.log
echo "Start:" $curDate > $log
# snakemake -s /work/sdukoldby/scripts/varscan2-joint-calling-wes-somatic-hg38.smk -j 999 --cluster "sbatch -A sdukoldby_slim --time 10:00:00" 2>&1 | tee -a $log
snakemake -s /work/sdukoldby/scripts/varscan2-joint-calling-deep-wes-somatic-hg38.smk -j 999 --cluster "sbatch -A sdukoldby_slim --time 10:00:00" 2>&1 | tee -a $log
echo "Slut:" `date +%F_%H-%M` >> $log
