#cd /scratch/sduhumac/kristina/data/genova/plasma/output/
cd /work/sduhumac/kristina/data/genova/tumor_blood/output/

sbatch /work/sduhumac/kristina/scripts/qualimap_fastqc.sh G35-1_illumina-truseq-genome_merged 
sbatch /work/sduhumac/kristina/scripts/qualimap_fastqc.sh G35-4_illumina-truseq-genome_merged


#for i in *_sorted_nodup.bam;
#do
#SAMPLEID=$(basename $i _sorted_nodup.bam)
#sbatch /work/sduhumac/kristina/scripts/qualimap_fastqc.sh $SAMPLEID
echo $SAMPLEID
#done
