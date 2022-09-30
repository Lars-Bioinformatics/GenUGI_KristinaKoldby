curDate=`date +"%F_%H-%M"`
log=log_mutect2-joint-calling-gatk4-wes-somatic-hg38-merged-samples_$curDate.log
echo "Start:" $curDate > $log
snakemake -s /work/sdukoldby/scripts/mutect2-joint-calling-gatk4-wes-somatic-hg38-merged-samples.smk -j 999 --cluster "sbatch -A sdukoldby_slim --time 10:00:00" 2>&1 | tee -a $log
echo "Slut:" `date +%F_%H-%M` >> $log
