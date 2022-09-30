
DAT="/work/sdukoldby/data/G45-2016_genugi/190325_A00653_0016_AHGL2LDSXX/BaseCalls/hg38/meerkat"

#SAMPLEID=$1
#SAMPLEID="G45-ECV2-4-blod_truseq-nano-genome_HGL2LDSXX_S26"
#SAMPLEID="G45-ECV2-8-blod_truseq-nano-genome_HGL2LDSXX_S29"
#SAMPLEID="G45-ECV2-29-blod_truseq-nano-genome_HGL2LDSXX_S30"
#SAMPLEID="G45-ECV2-31-blod_truseq-nano-genome_HGL2LDSXX_S27"
SAMPLEID="G45-ECV2-35-blod_truseq-nano-genome_HGL2LDSXX_S28"
#SAMPLEID="G45-ECV2-4-biopsi-H2_truseq-nano-genome_HGL2LDSXX_S24"
#SAMPLEID="G45-ECV2-8-biopsi-I2_truseq-nano-genome_HGL2LDSXX_S25"
#SAMPLEID="G45-ECV2-29-biopsi-C1_truseq-nano-genome_HGL2LDSXX_S21"
#SAMPLEID="G45-ECV2-31-biopsi-F1_truseq-nano-genome_HGL2LDSXX_S22"
#SAMPLEID="G45-ECV2-35-biopsi-G1_truseq-nano-genome_HGL2LDSXX_S23"

#SAMPLEID2="G45-ECV2-4-blod"
#SAMPLEID2="G45-ECV2-8-blod"
#SAMPLEID2="G45-ECV2-29-blod"
#SAMPLEID2="G45-ECV2-31-blod"
SAMPLEID2="G45-ECV2-35-blod"
#SAMPLEID2="G45-ECV2-4-biopsi-H2"
#SAMPLEID2="G45-ECV2-8-biopsi-I2"
#SAMPLEID2="G45-ECV2-29-biopsi-C1"
#SAMPLEID2="G45-ECV2-31-biopsi-F1"
#SAMPLEID2="G45-ECV2-35-biopsi-G1"

### Split -.variants files according to variant type ###

cd $DAT/$SAMPLEID
mkdir -p output

awk -v OFS='\t' '{if ($1 == "del") print $6,$7,$8,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11}' $SAMPLEID.recalibrated.variants > output/$SAMPLEID2.meerkat.DEL.txt
awk -v OFS='\t' '{if ($2<$3) print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14; else print $1,$3,$2,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14;}' output/$SAMPLEID2.meerkat.DEL.txt > output/$SAMPLEID2.meerkat.DEL2.txt
rm output/$SAMPLEID2.meerkat.DEL.txt

awk -v OFS='\t' '{if ($1 == "del_ins") print $6,$7,$8,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16}' $SAMPLEID.recalibrated.variants > output/$SAMPLEID2.meerkat.DEL_INS.txt
awk -v OFS='\t' '{if ($2<$3) print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19; else print $1,$3,$2,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19;}' output/$SAMPLEID2.meerkat.DEL_INS.txt > output/$SAMPLEID2.meerkat.DEL_INS2.txt
rm output/$SAMPLEID2.meerkat.DEL_INS.txt

awk -v OFS='\t' '{if ($1 == "del_insod") print $6,$7,$8,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16}' $SAMPLEID.recalibrated.variants > output/$SAMPLEID2.meerkat.DEL_INSOD.txt
awk -v OFS='\t' '{if ($2<$3) print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19; else print $1,$3,$2,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19;}' output/$SAMPLEID2.meerkat.DEL_INSOD.txt > output/$SAMPLEID2.meerkat.DEL_INSOD2.txt
rm output/$SAMPLEID2.meerkat.DEL_INSOD.txt

awk -v OFS='\t' '{if ($1 == "del_insou") print $6,$7,$8,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16}' $SAMPLEID.recalibrated.variants > output/$SAMPLEID2.meerkat.DEL_INSOU.txt
awk -v OFS='\t' '{if ($2<$3) print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12,$13,$14,$15,$16,$17,$18,$19; else print $1, $3, $2, $4, $5, $6, $7, $8, $9, $10, $11,$12,$13,$14,$15,$16,$17,$18,$19;}' output/$SAMPLEID2.meerkat.DEL_INSOU.txt > output/$SAMPLEID2.meerkat.DEL_INSOU2.txt
rm output/$SAMPLEID2.meerkat.DEL_INSOU.txt

awk -v OFS='\t' '{if ($1 == "del_inss") print $6,$7,$8,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16}' $SAMPLEID.recalibrated.variants > output/$SAMPLEID2.meerkat.DEL_INSS.txt
awk -v OFS='\t' '{if ($2<$3) print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19; else print $1,$3,$2,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19;}' output/$SAMPLEID2.meerkat.DEL_INSS.txt > output/$SAMPLEID2.meerkat.DEL_INSS2.txt
rm output/$SAMPLEID2.meerkat.DEL_INSS.txt

awk -v OFS='\t' '{if ($1 == "del_inssd") print $6,$7,$8,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16}' $SAMPLEID.recalibrated.variants > output/$SAMPLEID2.meerkat.DEL_INSSD.txt
awk -v OFS='\t' '{if ($2<$3) print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19; else print $1,$3,$2,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19;}' output/$SAMPLEID2.meerkat.DEL_INSSD.txt > output/$SAMPLEID2.meerkat.DEL_INSSD2.txt
rm output/$SAMPLEID2.meerkat.DEL_INSSD.txt

awk -v OFS='\t' '{if ($1 == "del_inssu") print $6,$7,$8,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16}' $SAMPLEID.recalibrated.variants > output/$SAMPLEID2.meerkat.DEL_INSSU.txt
awk -v OFS='\t' '{if ($2<$3) print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19; else print $1,$3,$2,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19;}' output/$SAMPLEID2.meerkat.DEL_INSSU.txt > output/$SAMPLEID2.meerkat.DEL_INSSU2.txt
rm output/$SAMPLEID2.meerkat.DEL_INSSU.txt

awk -v OFS='\t' '{if ($1 == "del_invers") print $6,$7,$8,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17}' $SAMPLEID.recalibrated.variants > output/$SAMPLEID2.meerkat.DEL_INVERS.txt
awk -v OFS='\t' '{if ($2<$3) print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20; else print $1,$3,$2,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20;}' output/$SAMPLEID2.meerkat.DEL_INVERS.txt > output/$SAMPLEID2.meerkat.DEL_INVERS2.txt
rm output/$SAMPLEID2.meerkat.DEL_INVERS.txt

awk -v OFS='\t' '{if ($1 == "inss") print $6,$7,$8,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16}' $SAMPLEID.recalibrated.variants > output/$SAMPLEID2.meerkat.INSS.txt
awk -v OFS='\t' '{if ($2<$3) print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19; else print $1,$3,$2,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19;}' output/$SAMPLEID2.meerkat.INSS.txt > output/$SAMPLEID2.meerkat.INSS2.txt
rm output/$SAMPLEID2.meerkat.INSS.txt

awk -v OFS='\t' '{if ($1 == "inssd") print $6,$7,$8,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16}' $SAMPLEID.recalibrated.variants > output/$SAMPLEID2.meerkat.INSSD.txt
awk -v OFS='\t' '{if ($2<$3) print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19; else print $1,$3,$2,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19;}' output/$SAMPLEID2.meerkat.INSSD.txt > output/$SAMPLEID2.meerkat.INSSD2.txt
rm output/$SAMPLEID2.meerkat.INSSD.txt

awk -v OFS='\t' '{if ($1 == "inssu") print $6,$7,$8,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12, $13,$14,$15,$16}' $SAMPLEID.recalibrated.variants > output/$SAMPLEID2.meerkat.INSSU.txt
awk -v OFS='\t' '{if ($2<$3) print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19; else print $1,$3,$2,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19;}' output/$SAMPLEID2.meerkat.INSSU.txt > output/$SAMPLEID2.meerkat.INSSU2.txt 
rm output/$SAMPLEID2.meerkat.INSSU.txt

awk -v OFS='\t' '{if ($1 == "insod") print $6,$7,$8,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16}' $SAMPLEID.recalibrated.variants > output/$SAMPLEID2.meerkat.INSOD.txt
awk -v OFS='\t' '{if ($2<$3) print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19; else print $1,$3,$2,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19;}' output/$SAMPLEID2.meerkat.INSOD.txt > output/$SAMPLEID2.meerkat.INSOD2.txt
rm output/$SAMPLEID2.meerkat.INSOD.txt

awk -v OFS='\t' '{if ($1 == "insou") print $6,$7,$8,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16}' $SAMPLEID.recalibrated.variants > output/$SAMPLEID2.meerkat.INSOU.txt
awk -v OFS='\t' '{if ($2<$3) print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19; else print $1,$3,$2,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19;}' output/$SAMPLEID2.meerkat.INSOU.txt > output/$SAMPLEID2.meerkat.INSOU2.txt
rm output/$SAMPLEID2.meerkat.INSOU.txt

awk -v OFS='\t' '{if ($1 == "invers") print $6,$7,$8,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' $SAMPLEID.recalibrated.variants > output/$SAMPLEID2.meerkat.INVERS.txt
awk -v OFS='\t' '{if ($2<$3) print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15; else print $1,$3,$2,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15;}' output/$SAMPLEID2.meerkat.INVERS.txt > output/$SAMPLEID2.meerkat.INVERS2.txt
rm output/$SAMPLEID2.meerkat.INVERS.txt

awk -v OFS='\t' '{if ($1 == "invers_f") print $6,$7,$8,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11}' $SAMPLEID.recalibrated.variants > output/$SAMPLEID2.meerkat.INVERS_f.txt
awk -v OFS='\t' '{if ($2<$3) print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14; else print $1,$3,$2,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14;}' output/$SAMPLEID2.meerkat.INVERS_f.txt > output/$SAMPLEID2.meerkat.INVERS_f2.txt
rm output/$SAMPLEID2.meerkat.INVERS_f.txt

awk -v OFS='\t' '{if ($1 == "invers_r") print $6,$7,$8,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11}' $SAMPLEID.recalibrated.variants > output/$SAMPLEID2.meerkat.INVERS_r.txt
awk -v OFS='\t' '{if ($2<$3) print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14; else print $1,$3,$2,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14;}' output/$SAMPLEID2.meerkat.INVERS_r.txt > output/$SAMPLEID2.meerkat.INVERS_r2.txt
rm output/$SAMPLEID2.meerkat.INVERS_r.txt 

awk -v OFS='\t' '{if ($1 == "tandem_dup") print $6,$7,$8,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11}' $SAMPLEID.recalibrated.variants > output/$SAMPLEID2.meerkat.TANDEM_DUP.txt
awk -v OFS='\t' '{if ($2<$3) print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14; else print $1,$3,$2,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14;}' output/$SAMPLEID2.meerkat.TANDEM_DUP.txt > output/$SAMPLEID2.meerkat.TANDEM_DUP2.txt
rm output/$SAMPLEID2.meerkat.TANDEM_DUP.txt

awk -v OFS='\t' '{if ($1 == "transl_inter") print $6,$7,$8,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13}' $SAMPLEID.recalibrated.variants > output/$SAMPLEID2.meerkat.TRANSL_INTER.txt
awk -v OFS='\t' '{if ($2<$3) print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16; else print $1,$3,$2,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16;}' output/$SAMPLEID2.meerkat.TRANSL_INTER.txt > output/$SAMPLEID2.meerkat.TRANSL_INTER2.txt
rm output/$SAMPLEID2.meerkat.TRANSL_INTER.txt


