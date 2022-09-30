MER="/data/mark/tools/Meerkat2"
ALIGNED="/data/mark/tools/breakpoint_aln"
DAT="/data/mark/tools/Meerkat2/data"
BWA="/data/mark/tools/bwa062"
SAM="/data/mark/tools/samtools"
BLAST="/data/mark/tools/blast-2.2.24/bin"
REF="/data/mark/tools/ref"

ALIGNED="/data/mark/tools/breakpoint_aln"

for i in *.bam;
do

newfile=$(basename $i .bam)
echo "$newfile"

cp $ALIGNED/${newfile}.bam $DAT/
cp $ALIGNED/${newfile}.bai $DAT/

cd $DAT

perl $MER/scripts/pre_process.pl -b ${newfile}.bam -I $REF/hg19.fa -A $REF/hg19.fa.fai -W $BWA -S $SAM -t 24

perl $MER/scripts/meerkat.pl -b ${newfile}.bam -F $RAP/hg19.fa -W $BWA -B $BLAST -S $SAM -t 24

perl $MER/scripts/mechanism.pl -b ${newfile}.bam -R $REF/rmsk-hg19.txt

awk -v OFS='\t' '{if ($1 == "del") print $6,$7,$8,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11}' $DAT/sample.sorted.nodup.variants > $DAT/smp_DEL.txt
awk -v OFS='\t' '{if ($2<$3) print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14; else print $1,$3,$2,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14;}' $DAT/smp_DEL.txt > $DAT/smp_DEL2.txt
rm $DAT/smp_DEL.txt
awk -v OFS='\t' '{if ($1 == "del_ins") print $6,$7,$8,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16}' $DAT/sample.sorted.nodup.variants > $DAT/smp_DEL_INS.txt
awk -v OFS='\t' '{if ($2<$3) print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19; else print $1,$3,$2,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19;}' $DAT/smp_DEL_INS.txt > $DAT/smp_DEL_INS2.txt
rm $DAT/smp_DEL_INS.txt
awk -v OFS='\t' '{if ($1 == "del_insod") print $6,$7,$8,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16}' $DAT/sample.sorted.nodup.variants > $DAT/smp_DEL_INSOD.txt
awk -v OFS='\t' '{if ($2<$3) print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19; else print $1,$3,$2,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19;}' $DAT/smp_DEL_INSOD.txt > $DAT/smp_DEL_INSOD2.txt
rm $DAT/smp_DEL_INSOD.txt
awk -v OFS='\t' '{if ($1 == "del_insou") print $6,$7,$8,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16}' $DAT/sample.sorted.nodup.variants > $DAT/smp_DEL_INSOU.txt
awk -v OFS='\t' '{if ($2<$3) print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12,$13,$14,$15,$16,$17,$18,$19; else print $1, $3, $2, $4, $5, $6, $7, $8, $9, $10, $11,$12,$13,$14,$15,$16,$17,$18,$19;}' $DAT/smp_DEL_INSOU.txt > $DAT/smp_DEL_INSOU2.txt
rm $DAT/smp_DEL_INSOU.txt
awk -v OFS='\t' '{if ($1 == "del_inss") print $6,$7,$8,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16}' $DAT/sample.sorted.nodup.variants > $DAT/smp_DEL_INSS.txt
awk -v OFS='\t' '{if ($2<$3) print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19; else print $1,$3,$2,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19;}' $DAT/smp_DEL_INSS.txt > $DAT/smp_DEL_INSS2.txt
rm $DAT/smp_DEL_INSS.txt
awk -v OFS='\t' '{if ($1 == "del_inssd") print $6,$7,$8,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16}' $DAT/sample.sorted.nodup.variants > $DAT/smp_DEL_INSSD.txt
awk -v OFS='\t' '{if ($2<$3) print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19; else print $1,$3,$2,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19;}' $DAT/smp_DEL_INSSD.txt > $DAT/smp_DEL_INSSD2.txt
rm $DAT/smp_DEL_INSSD.txt
awk -v OFS='\t' '{if ($1 == "del_inssu") print $6,$7,$8,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16}' $DAT/sample.sorted.nodup.variants > $DAT/smp_DEL_INSSU.txt
awk -v OFS='\t' '{if ($2<$3) print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19; else print $1,$3,$2,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19;}' $DAT/smp_DEL_INSSU.txt > $DAT/smp_DEL_INSSU2.txt
rm $DAT/smp_DEL_INSSU.txt
awk -v OFS='\t' '{if ($1 == "del_invers") print $6,$7,$8,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17}' $DAT/sample.sorted.nodup.variants > $DAT/smp_DEL_INVERS.txt
awk -v OFS='\t' '{if ($2<$3) print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20; else print $1,$3,$2,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20;}' $DAT/smp_DEL_INVERS.txt > $DAT/smp_DEL_INVERS2.txt
rm $DAT/smp_DEL_INVERS.txt
awk -v OFS='\t' '{if ($1 == "inss") print $6,$7,$8,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16}' $DAT/sample.sorted.nodup.variants > $DAT/smp_INSS.txt
awk -v OFS='\t' '{if ($2<$3) print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19; else print $1,$3,$2,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19;}' $DAT/smp_INSS.txt > $DAT/smp_INSS2.txt
rm $DAT/smp_INSS.txt
awk -v OFS='\t' '{if ($1 == "inssd") print $6,$7,$8,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16}' $DAT/sample.sorted.nodup.variants > $DAT/smp_INSSD.txt
awk -v OFS='\t' '{if ($2<$3) print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19; else print $1,$3,$2,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19;}' $DAT/smp_INSSD.txt > $DAT/smp_INSSD2.txt
rm $DAT/smp_INSSD.txt
awk -v OFS='\t' '{if ($1 == "inssu") print $6,$7,$8,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12, $13,$14,$15,$16}' $DAT/sample.sorted.nodup.variants > $DAT/smp_INSSU.txt
awk -v OFS='\t' '{if ($2<$3) print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19; else print $1,$3,$2,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19;}' $DAT/smp_INSSU.txt > $DAT/smp_INSSU2.txt 
awk -v OFS='\t' '{if ($1 == "insod") print $6,$7,$8,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16}' $DAT/sample.sorted.nodup.variants > $DAT/smp_INSOD.txt

awk -v OFS='\t' '{if ($1 == "insod") print $6,$7,$8,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16}' $DAT/sample.sorted.nodup.variants > $DAT/smp_INSOD.txt

awk -v OFS='\t' '{if ($2<$3) print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19; else print $1,$3,$2,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19;}' $DAT/smp_INSOD.txt > $DAT/smp_INSOD2.txt

awk -v OFS='\t' '{if ($1 == "insou") print $6,$7,$8,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16}' $DAT/sample.sorted.nodup.variants > $DAT/smp_INSOU.txt
awk -v OFS='\t' '{if ($2<$3) print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19; else print $1,$3,$2,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19;}' $DAT/smp_INSOU.txt > $DAT/smp_INSOU2.txt
rm $DAT/smp_INSOU.txt
awk -v OFS='\t' '{if ($1 == "invers") print $6,$7,$8,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' $DAT/sample.sorted.nodup.variants > $DAT/smp_INVERS.txt
awk -v OFS='\t' '{if ($2<$3) print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15; else print $1,$3,$2,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15;}' $DAT/smp_INVERS.txt > $DAT/smp_INVERS2.txt
rm $DAT/smp_INVERS.txt
awk -v OFS='\t' '{if ($1 == "invers_f") print $6,$7,$8,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11}' $DAT/sample.sorted.nodup.variants > $DAT/smp_INVERS_f.txt
awk -v OFS='\t' '{if ($2<$3) print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14; else print $1,$3,$2,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14;}' $DAT/smp_INVERS_f.txt > $DAT/smp_INVERS_f2.txt
rm $DAT/smp_INVERS_f.txt
awk -v OFS='\t' '{if ($1 == "invers_r") print $6,$7,$8,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11}' $DAT/sample.sorted.nodup.variants > $DAT/smp_INVERS_r.txt
awk -v OFS='\t' '{if ($2<$3) print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14; else print $1,$3,$2,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14;}' $DAT/smp_INVERS_r.txt > $DAT/smp_INVERS_r2.txt
rm $DAT/smp_INVERS_r.txt 
awk -v OFS='\t' '{if ($1 == "tandem_dup") print $6,$7,$8,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11}' $DAT/sample.sorted.nodup.variants > $DAT/smp_TANDEM_DUP.txt
awk -v OFS='\t' '{if ($2<$3) print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14; else print $1,$3,$2,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14;}' $DAT/smp_TANDEM_DUP.txt > $DAT/smp_TANDEM_DUP2.txt
rm $DAT/smp_TANDEM_DUP.txt
awk -v OFS='\t' '{if ($1 == "transl_inter") print $6,$7,$8,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13}' $DAT/sample.sorted.nodup.variants > $DAT/smp_TRANSL_INTER.txt
awk -v OFS='\t' '{if ($2<$3) print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16; else print $1,$3,$2,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16;}' $DAT/smp_TRANSL_INTER.txt > $DAT/smp_TRANSL_INTER2.txt

rm $DAT/smp_TRANSL_INTER.txt

done



perl $MER/scripts/somatic_sv.pl -i $RAP/sample.sorted.nodup.variants -o $RAP/somatic_filtered.sample.sorted.nodup.variants

perl $MER/scripts/meerkat2vcf.pl -i $RAP/somatic_filtered.sample.sorted.nodup.variants -H /data/mark/tools/Meerkat/Meerkat.example/headerfile -F /data/mark/tools/ref/hg19.fasta -o $RAP/somatic_filtered.sample.sorted.nodup.variants.vcf

perl $MER/scripts/fusions.pl -i $RAP/somatic_filtered.sample.sorted.nodup.variant -G /data/mark/tools/RAPTR-SV/genes.bed
