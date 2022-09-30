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
--readFilesIn /gpfs/gss1/work/sduhumac/data/jeanette/G56-sampleA1_rnaseq_HGTHHDSXX_R1.fastq.gz /gpfs/gss1/work/sduhumac/data/jeanette/G56-sampleA1_rnaseq_HGTHHDSXX_R2.fastq.gz \
--outFileNamePrefix A1 \
--quantMode GeneCounts \
--readFilesCommand zcat \
--chimSegmentMin 20 \
--chimOutType Junctions \
--outSAMtype BAM SortedByCoordinate

cd /gpfs/gss1/work/sduhumac/data/jeanette
/gpfs/gss1/work/sduhumac/tools/miniconda3/bin/STAR \
--runThreadN 24 \
--genomeDir /gpfs/gss1/work/sduhumac/tools/cellranger-2.1.1/refdata-cellranger-GRCh38-3.0.0/star/ \
--sjdbGTFfile /gpfs/gss1/work/sduhumac/tools/cellranger-2.1.1/refdata-cellranger-GRCh38-3.0.0/genes/genes.gtf \
--readFilesIn /gpfs/gss1/work/sduhumac/data/jeanette/G56-sampleA8_rnaseq_HGTHHDSXX_R1.fastq.gz /gpfs/gss1/work/sduhumac/data/jeanette/G56-sampleA8_rnaseq_HGTHHDSXX_R2.fastq.gz \
--outFileNamePrefix A8 \
--quantMode GeneCounts \
--readFilesCommand zcat \
--chimSegmentMin 20 \
--chimOutType Junctions \
--outSAMtype BAM SortedByCoordinate

cd /gpfs/gss1/work/sduhumac/data/jeanette
/gpfs/gss1/work/sduhumac/tools/miniconda3/bin/STAR \
--runThreadN 24 \
--genomeDir /gpfs/gss1/work/sduhumac/tools/cellranger-2.1.1/refdata-cellranger-GRCh38-3.0.0/star/ \
--sjdbGTFfile /gpfs/gss1/work/sduhumac/tools/cellranger-2.1.1/refdata-cellranger-GRCh38-3.0.0/genes/genes.gtf \
--readFilesIn /gpfs/gss1/work/sduhumac/data/jeanette/G56-sampleA9_rnaseq_HGTHHDSXX_R1.fastq.gz /gpfs/gss1/work/sduhumac/data/jeanette/G56-sampleA9_rnaseq_HGTHHDSXX_R2.fastq.gz \
--outFileNamePrefix A9 \
--quantMode GeneCounts \
--readFilesCommand zcat \
--chimSegmentMin 20 \
--chimOutType Junctions \
--outSAMtype BAM SortedByCoordinate

cd /gpfs/gss1/work/sduhumac/data/jeanette
/gpfs/gss1/work/sduhumac/tools/miniconda3/bin/STAR \
--runThreadN 24 \
--genomeDir /gpfs/gss1/work/sduhumac/tools/cellranger-2.1.1/refdata-cellranger-GRCh38-3.0.0/star/ \
--sjdbGTFfile /gpfs/gss1/work/sduhumac/tools/cellranger-2.1.1/refdata-cellranger-GRCh38-3.0.0/genes/genes.gtf \
--readFilesIn /gpfs/gss1/work/sduhumac/data/jeanette/G56-sampleB1_rnaseq_HGTHHDSXX_R1.fastq.gz /gpfs/gss1/work/sduhumac/data/jeanette/G56-sampleB1_rnaseq_HGTHHDSXX_R2.fastq.gz \
--outFileNamePrefix B1 \
--quantMode GeneCounts \
--readFilesCommand zcat \
--chimSegmentMin 20 \
--chimOutType Junctions \
--outSAMtype BAM SortedByCoordinate

cd /gpfs/gss1/work/sduhumac/data/jeanette
/gpfs/gss1/work/sduhumac/tools/miniconda3/bin/STAR \
--runThreadN 24 \
--genomeDir /gpfs/gss1/work/sduhumac/tools/cellranger-2.1.1/refdata-cellranger-GRCh38-3.0.0/star/ \
--sjdbGTFfile /gpfs/gss1/work/sduhumac/tools/cellranger-2.1.1/refdata-cellranger-GRCh38-3.0.0/genes/genes.gtf \
--readFilesIn /gpfs/gss1/work/sduhumac/data/jeanette/G56-sampleB2_rnaseq_HGTHHDSXX_R1.fastq.gz /gpfs/gss1/work/sduhumac/data/jeanette/G56-sampleB2_rnaseq_HGTHHDSXX_R2.fastq.gz \
--outFileNamePrefix B2 \
--quantMode GeneCounts \
--readFilesCommand zcat \
--chimSegmentMin 20 \
--chimOutType Junctions \
--outSAMtype BAM SortedByCoordinate


