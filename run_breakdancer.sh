cd /work/sdukoldby/data/G45-2016_genugi/190325_A00653_0016_AHGL2LDSXX/BaseCalls/hg38

for i in *.recalibrated.bam;
do
SAMPLEID=$(basename $i .recalibrated.bam)
echo $SAMPLEID
sbatch /work/sdukoldby/scripts/breakdancer.sh $SAMPLEID
done
