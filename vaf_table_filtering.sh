
FILTEREDVCF="ECV2-4_tagseq-medexome_somatic_mutect2_vaf_unfiltered.txt"
PASSFILTER="ECV2-4_tagseq-medexome_somatic_mutect2_vaf_filtered_pass.txt"
ADFILTER="ECV2-4_tagseq-medexome_somatic_mutect2_vaf_filtered_ADmin3.txt"
minAD=3


### Filtering: FILTER = PASS, keep header ###

# awk '{if(NR==1 || $3=="PASS") { print $0 }}' $FILTEREDVCF > $PASSFILTER


### Filtering on AD ###

# awk '{if($4>=$minAD || \
# $5>=$minAD || \
# $7>=$minAD || \
# $8>=$minAD || \
# $9>=$minAD || \
# $10>=$minAD || \
# $11>=$minAD || \
# $12>=$minAD) \
# {print $0}}' \
# $PASSFILTER
# > $ADFILTER

# echo bcftools filter --include 'FILTER=="PASS" && FORMAT/AD[*:1]>=3' ECV2-4_tagseq-medexome_somatic_mutect2_filtered.vcf
# echo ECV2-4_tagseq-medexome_somatic_mutect2_vaf_filtered_pass_ADmin${minAD}.txt

bcftools filter --include 'FILTER=="PASS" && FORMAT/AD[*:1]>=3' ECV2-4_tagseq-medexome_somatic_mutect2_filtered.vcf |\
bcftools query -H -f '%CHROM\t%POS\t%FILTER\t[%AF\t][%AD{1}\t][%DP\t]\n' > ECV2-4_tagseq-medexome_somatic_mutect2_vaf_filtered_pass_ADmin3.txt

###  Fjern linjer hvor dybde i blodprøve (nr 2 i 0-baseret tælling) er mindre end 20.
# && FORMAT/DP[2] >= 20
