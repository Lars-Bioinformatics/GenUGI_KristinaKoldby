f_cur=""
f_cur_path=""
for i in *normal*R1*gz; do
  f="$( cut -d '_' -f 1 <<< "$i" )"
  echo $f
  if [[ $f != ${f_cur} ]]; then
    if [[ $f_cur == "" ]]; then
      f_cur=$f
      f_cur_path=$i
    else
      cp ${f_cur_path} exome_fastq_merged/${f}_normal_tagseq-medexome_R1_001.fastq.gz
      f_cur=$f
      f_cur_path=$i
    fi
  else
    zcat ${f_cur_path} $i | gzip -c > exome_fastq_merged/${f}_normal_tagseq-medexome_R1_001.fastq.gz
    f_cur=""
  fi
done &

f_cur2=""
f_cur2_path=""
for i in *normal*R2*gz; do
  f="$( cut -d '_' -f 1 <<< "$i" )"
  echo $f
if [[ $f != ${f_cur2} ]]; then
  if [[ $f_cur2 == "" ]]; then
      f_cur2=$f
      f_cur2_path=$i
    else
      cp ${f_cur2_path} exome_fastq_merged/${f}_normal_tagseq-medexome_R2_001.fastq.gz
      f_cur2=$f
      f_cur2_path=$i
    fi
  else
    zcat ${f_cur2_path} $i | gzip -c > exome_fastq_merged/${f}_normal_tagseq-medexome_R2_001.fastq.gz
    f_cur2=""
  fi
done &

f_cur3=""
f_cur3_path=""
for i in *tumor*R1*gz; do
  f="$( cut -d '_' -f 1 <<< "$i" )"
  echo $f
if [[ $f != ${f_cur3} ]]; then
  if [[ $f_cur3 == "" ]]; then
      f_cur3=$f
      f_cur3_path=$i
    else
      cp ${f_cur3_path} exome_fastq_merged/${f}_tumor_tagseq-medexome_R1_001.fastq.gz
      f_cur3=$f
      f_cur3_path=$i
    fi
  else
    zcat ${f_cur3_path} $i | gzip -c > exome_fastq_merged/${f}_tumor_tagseq-medexome_R1_001.fastq.gz
    f_cur3=""
  fi
done &

f_cur4=""
f_cur4_path=""
for i in *tumor*R2*gz; do
  f="$( cut -d '_' -f 1 <<< "$i" )"
  echo $f
if [[ $f != ${f_cur4} ]]; then
  if [[ $f_cur4 == "" ]]; then
      f_cur4=$f
      f_cur4_path=$i
    else
      cp ${f_cur4_path} exome_fastq_merged/${f}_tumor_tagseq-medexome_R2_001.fastq.gz
      f_cur4=$f
      f_cur4_path=$i
    fi
  else
    zcat ${f_cur4_path} $i | gzip -c > exome_fastq_merged/${f}_tumor_tagseq-medexome_R2_001.fastq.gz
    f_cur4=""
  fi
done


# for i in *normal*R2*gz; do
#   f="$( cut -d '_' -f 1 <<< "$i" )"
#   echo $f
#   less $i >> exome_fastq_merged/${f}_normal_tagseq-medexome_R2_001.fastq
# done &
#
# for i in *tumor*R1*gz; do
#   f="$( cut -d '_' -f 1 <<< "$i" )"
#   echo $f
#   less $i >> exome_fastq_merged/${f}_tumor_tagseq-medexome_R1_001.fastq
# done &
#
# for i in *tumor*R2*gz; do
#   f="$( cut -d '_' -f 1 <<< "$i" )"
#   echo $f
#   less $i >> exome_fastq_merged/${f}_tumor_tagseq-medexome_R2_001.fastq
# done
