curDate=`date +"%F_%H-%M"`
log=log_mutect2-gatk4-wes-somatic-hg38_$curDate.log
echo "Start:" $curDate > $log
snakemake -s /work/sdukoldby/scripts/mutect2-gatk4-wes-somatic-hg38.smk -j 999 --cluster "sbatch -A sdukoldby_slim --time 10:00:00" 2>&1 | tee -a $log
echo "Slut:" `date +%F_%H-%M` >> $log
