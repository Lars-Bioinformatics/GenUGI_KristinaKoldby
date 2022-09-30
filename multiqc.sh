#!/bin/sh
#
#SBATCH --account sduhumac_slim      # account
#SBATCH --nodes 1                 # number of nodes
#SBATCH --time 24:00:00            # max time (HH:MM:SS)
#SBATCH --output=/work/sduhumac/kristina/scripts/slurm_log/multiqc_%j.out

set -x
echo $PWD

#cd /work/sduhumac/kristina/data/genova/tumor_blood/
cd /scratch/sduhumac/kristina/data/genova/plasma/output/qc/

### Collect QC reports for all samples with multiqc
multiqc . /work/sduhumac/kristina/data/genova/tumor_blood/output/qc/ \
--config /work/sduhumac/kristina/scripts/multiqc_config.yaml \
--outdir /scratch/sduhumac/kristina/data/genova/plasma/output/qc/