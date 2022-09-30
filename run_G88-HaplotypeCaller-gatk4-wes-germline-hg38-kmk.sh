curDate=`date +"%F_%H-%M"`
log=log_HaplotypeCaller-gatk4-wes-germline-hg38_$curDate.log
echo "Start:" $curDate > $log
snakemake -s /work/sdukoldby/scripts/G88-HaplotypeCaller-gatk4-wes-germline-hg38-kmk.smk -j 999 --cluster "sbatch -A sdukoldby_slim --time 24:00:00" 2>&1 | tee -a $log
echo "Slut:" `date +%F_%H-%M` >> $log
