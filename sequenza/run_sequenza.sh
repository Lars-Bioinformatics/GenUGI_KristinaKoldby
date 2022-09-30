# snakemake -s /work/sdukoldby/scripts/sequenza/sequenza.smk -j 999 --cluster "sbatch -A sdukoldby_slim --time 24:00:00"
# snakemake -s /work/sdukoldby/scripts/sequenza/sequenza_deep-seq.smk -j 999 --cluster "sbatch -A sdukoldby_fat --time 01:00:00"
# snakemake -s /work/sdukoldby/scripts/sequenza/sequenza_wgs.smk -j 999 --cluster "sbatch --qos=long -A sdukoldby_fat --time 120:00:00"
snakemake -s /work/sdukoldby/scripts/sequenza/sequenza_wgs.smk -j 999 --cluster "sbatch -A sdukoldby_slim --time 24:00:00"
