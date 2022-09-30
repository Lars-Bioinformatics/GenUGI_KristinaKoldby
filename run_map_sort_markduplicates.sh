#cd /work/sduhumac/kristina/data/genova/tumor_blood/
cd /work/sduhumac/kristina/data/genova/plasma/

for i in *_R1.fastq.gz;
do
SAMPLEID=$(basename $i _R1.fastq.gz)
#echo $SAMPLEID
if [[ $SAMPLEID == "G35-2"* ]]; then
#LOGFILE="/scratch/sduhumac/kristina/data/genova/plasma/output/logfiles/$SAMPLEID_log_merge_bam_markduplicates.txt"
sbatch /work/sduhumac/kristina/scripts/map_sort_markduplicates.sh $SAMPLEID
echo $SAMPLEID
fi
done
