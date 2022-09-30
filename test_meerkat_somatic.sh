SAMPLEID="G45-ECV2-4-biopsi-H2_truseq-nano-genome_HGL2LDSXX_S24"
MER="/work/sdukoldby/tools/Meerkat"
SCRIPTS="/work/sdukoldby/tools/Meerkat/scripts"
DAT="/work/sdukoldby/data/G45-2016_genugi/190325_A00653_0016_AHGL2LDSXX/BaseCalls/hg38/meerkat"
OUT="/work/sdukoldby/data/G45-2016_genugi/190325_A00653_0016_AHGL2LDSXX/BaseCalls/hg38/meerkat"
BLAST="/work/sdukoldby/tools/ncbi-blast-2.9.0+/bin"
REF="/work/sdukoldby/resources/hg38"
BIODB="/work/sdukoldby/tools/miniconda3/lib/site_perl/5.26.2/Bio/DB"

perl $MER/scripts/somatic_sv.pl -i $DAT/$SAMPLEID/$SAMPLEID.recalibrated.variants -o $DAT/$SAMPLEID/$SAMPLEID.somatic_filtered.variants -F $DAT/normal_discord_files/ -R $REF/rmsk-hg38.txt

