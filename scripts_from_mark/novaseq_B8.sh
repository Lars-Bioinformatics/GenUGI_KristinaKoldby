#!/bin/sh
#
#SBATCH --account sdukoldby_fat      # account
#SBATCH --nodes 1                # number of nodes
#SBATCH --time 24:00:00            # max time (HH:MM:SS)

cd /gpfs/gss1/work/sduhumac/data/jeanette
/gpfs/gss1/work/sduhumac/tools/miniconda3/bin/STAR \
--runThreadN 24 \
--genomeDir /gpfs/gss1/work/sduhumac/tools/cellranger-2.1.1/refdata-cellranger-GRCh38-3.0.0/star/ \
--sjdbGTFfile /gpfs/gss1/work/sduhumac/tools/cellranger-2.1.1/refdata-cellranger-GRCh38-3.0.0/genes/genes.gtf \
--readFilesIn /gpfs/gss1/work/sduhumac/data/jeanette/G56-sampleB8_rnaseq_HGTHHDSXX_R1.fastq.gz /gpfs/gss1/work/sduhumac/data/jeanette/G56-sampleB8_rnaseq_HGTHHDSXX_R2.fastq.gz \
--outFileNamePrefix B8 \
--quantMode GeneCounts \
--readFilesCommand zcat \
--chimSegmentMin 20 \
--chimOutType Junctions \
--outReadsUnmapped None
--twopassMode Basic
--outSAMstrandField intronMotif \
--outSAMtype BAM SortedByCoordinate\
--outSAMunmapped Within \
--chimSegmentMin 12 \  # ** essential to invoke chimeric read detection & reporting **
--chimJunctionOverhangMin 12 \
--chimOutJunctionFormat 1 \   # **essential** includes required metadata in Chimeric.junction.out file.
--alignSJDBoverhangMin 10 \
--alignMatesGapMax 100000 \   # avoid readthru fusions within 100k
--alignIntronMax 100000 \
--alignSJstitchMismatchNmax 5 -1 5 5 \   # settings improved certain chimera detections
--outSAMattrRGline ID:GRPundef \
--chimMultimapScoreRange 3 \
--chimScoreJunctionNonGTAG -4 \
--chimMultimapNmax 20 \
--chimNonchimScoreDropMin 10 \
--peOverlapNbasesMin 12 \
--peOverlapMMp 0.1 

