cd /work/sdukoldby/data/testdata_meerkat/data

for i in *.recalibrated.bam;
do
SAMPLEID=$(basename $i .recalibrated.bam)
echo $SAMPLEID
sbatch /work/sdukoldby/scripts/Meerkat.sh $SAMPLEID
done
