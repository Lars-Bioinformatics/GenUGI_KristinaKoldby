WORK=/work/sdukoldby/data/G45-2016_genugi/wgs_bam
REF=${WORK}/reference_files_GRCh38/archives

python2 ~/udocker-1.1.4/udocker run --user=pi --volume=$WORK:$WORK cgpwgs \
    ds-cgpwgs.pl \
    -r $REF/core_ref_GRCh38_hla_decoy_ebv.tar.gz \
    -a $REF/VAGrENT_ref_GRCh38_hla_decoy_ebv_ensembl_91.tar.gz \
    -si $REF/SNV_INDEL_ref_GRCh38_hla_decoy_ebv-fragment.tar.gz \
    -cs $REF/CNV_SV_ref_GRCh38_hla_decoy_ebv_brass6+.tar.gz \
    -qc $REF/qcGenotype_GRCh38_hla_decoy_ebv.tar.gz \
    -t $WORK/cgp_out2/old_cgpmap/G45-ECV2-31-biopsi-F1_truseq-nano-genome_HGL2LDSXX_S22/G45-ECV2-31-biopsi-F1_truseq-nano-genome_HGL2LDSXX_S22.bam \
    -tidx $WORK/cgp_out2/old_cgpmap/G45-ECV2-31-biopsi-F1_truseq-nano-genome_HGL2LDSXX_S22/G45-ECV2-31-biopsi-F1_truseq-nano-genome_HGL2LDSXX_S22.bam.bai \
    -n $WORK/cgp_out2/old_cgpmap/G45-ECV2-31-blod_truseq-nano-genome_HGL2LDSXX_S27/G45-ECV2-31-blod_truseq-nano-genome_HGL2LDSXX_S27.bam \
    -nidx $WORK/cgp_out2/old_cgpmap/G45-ECV2-31-blod_truseq-nano-genome_HGL2LDSXX_S27/G45-ECV2-31-blod_truseq-nano-genome_HGL2LDSXX_S27.bam.bai \
    -e HLA* \
    -o $WORK/test_cgpwgs
    