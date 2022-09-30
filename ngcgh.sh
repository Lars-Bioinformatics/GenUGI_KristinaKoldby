#!/bin/sh
#
#SBATCH --account sdukoldby_slim      # account
#SBATCH --nodes 1                 # number of nodes
#SBATCH --time 24:00:00            # max time (HH:MM:SS)

conda activate py2

## Form:
##ngCGH -w WINDOWSIZE (default 1000) -o OUTPUT normalbam tumorbam

#ngCGH -w 10000 -o ngcgh/G45-ECV2-4-biopsi-H2_truseq-nano-genome.ngcgh G45-ECV2-4-blod_truseq-nano-genome_HGL2LDSXX_S26.recalibrated.bam G45-ECV2-4-biopsi-H2_truseq-nano-genome_HGL2LDSXX_S24.recalibrated.bam -t 24

#ngCGH -w 10000 -o ngcgh/G45-ECV2-8-biopsi-I2_truseq-nano-genome.ngcgh G45-ECV2-8-blod_truseq-nano-genome_HGL2LDSXX_S29.recalibrated.bam G45-ECV2-8-biopsi-I2_truseq-nano-genome_HGL2LDSXX_S25.recalibrated.bam -t 24

#ngCGH -w 10000 -o ngcgh/G45-ECV2-29-biopsi-C1_truseq-nano-genome.ngcgh G45-ECV2-29-blod_truseq-nano-genome_HGL2LDSXX_S30.recalibrated.bam G45-ECV2-29-biopsi-C1_truseq-nano-genome_HGL2LDSXX_S21.recalibrated.bam -t 24

#ngCGH -w 10000 -o ngcgh/G45-ECV2-31-biopsi-F1_truseq-nano-genome.ngcgh G45-ECV2-31-blod_truseq-nano-genome_HGL2LDSXX_S27.recalibrated.bam G45-ECV2-31-biopsi-F1_truseq-nano-genome_HGL2LDSXX_S22.recalibrated.bam -t 24

ngCGH -w 10000 -o ngcgh/G45-ECV2-35-biopsi-G1_truseq-nano-genome.ngcgh G45-ECV2-35-blod_truseq-nano-genome_HGL2LDSXX_S28.recalibrated.bam G45-ECV2-35-biopsi-G1_truseq-nano-genome_HGL2LDSXX_S23.recalibrated.bam -t 24

        
  
       
   