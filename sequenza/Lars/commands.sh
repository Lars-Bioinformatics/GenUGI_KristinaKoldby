#! /bin/bash
#
#SBATCH --account sduvarcall_slim		# account
#SBATCH --nodes 1						# number of nodes
#SBATCH --time 24:00:00					# max time (HH:MM:SS)

source activate py2

normal=$1
tumor=$2

nname=`basename ${normal%_nimblegen*}`
tname=`basename ${tumor%_nimblegen*}`

# echo $normal $tumor
echo $nname ${tname}

### Already in path
# SEQUENZA_UTILS=~/miniconda3/lib/R/library/sequenza/exec/sequenza-utils.py

out_dir=Sequenza_${nname}_${tname}/
mkdir -p $out_dir

seqz_data=${out_dir}${nname}_${tname}.seqz.gz
seqz_data_bin500=${out_dir}${nname}_${tname}.seqz.bin500.gz

REF=/work/sduvarcall/bwa-0.7.13/reference

# Generate a genome-wide GC content file (only needs to be done once for each reference)
## sequenza−utils.py GC−windows −w 50 human_g1k_v37_decoy.fasta | gzip > Homo_sapiens.GRCh37.gc50Base.txt.gz

# Which build to use? Mostly GRCh37
##-gc ${REF}/hg19.gc5Base.txt.gz \
##-gc ${REF}/Homo_sapiens.GRCh37.gc50Base.txt.gz \

sequenza-utils bam2seqz \
	-n $normal \
	-t $tumor \
	--fasta ${REF}/human_g1k_v37_decoy.fasta \
	-gc ${REF}/Homo_sapiens.GRCh37.gc50Base.txt.gz \
	-q 20 \
	-N 20 | gzip > $seqz_data
echo 'Done with bam2seqz'


sequenza-utils seqz_binning \
	-w 500 \
	-s $seqz_data \
	| bgzip -c > $seqz_data_bin500
echo 'Done with seqz-binning'

Rscript --vanilla sequenza.R $seqz_data_bin500 ${nname}_${tname} $out_dir
# Rscript --vanilla sequenza.R "${nname}_${tname}.seqz.bin500.noGL.noMT.gz" ${nname}_${tname} $out_dir

sed -i 's/"//g' ${out_dir}Sequenza_${nname}_${tname}_segments.txt
sed -i 's/"//g' ${out_dir}Sequenza_${nname}_${tname}_mutations.txt

