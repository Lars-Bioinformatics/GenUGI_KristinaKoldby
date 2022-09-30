#!/bin/sh
#
#SBATCH --account sdukoldby_fat      # account
#SBATCH --nodes 1                 # number of nodes
#SBATCH --time 02:00:00            # max time (HH:MM:SS)

# bcftools reheader -f /work/sdukoldby/resources/hg38/Homo_sapiens_assembly38.fasta.fai -o ECV2-4_cns_varscan2_somatic_pon-filtered_strict_newheader.vcf ECV2-4_cns_varscan2_somatic_pon-filtered_strict.vcf

gatk Funcotator --variant ECV2-4_cns_varscan2_somatic_pon-filtered_strict_newheader.vcf --reference /work/sdukoldby/resources/hg38/Homo_sapiens_assembly38.fasta --ref-version hg38 --data-sources-path /work/sdukoldby/resources/funcotator_dataSources.v1.6.20190124s --output ECV2-4_cns_varscan2_somatic_pon-filtered_strict_funcotated.vcf --output-file-format VCF
# gatk Funcotator --variant ECV2-4_cns_varscan2_somatic_pon-filtered_strict.vcf --disable-sequence-dictionary-validation --reference /work/sdukoldby/resources/hg38/Homo_sapiens_assembly38.fasta --ref-version hg38 --data-sources-path /work/sdukoldby/resources/funcotator_dataSources.v1.6.20190124s --output ECV2-4_cns_varscan2_somatic_pon-filtered_strict_funcotated.vcf --output-file-format VCF
# gatk Funcotator --variant ECV2-4_somatic_mutect2_filtered.vcf --reference /work/sdukoldby/resources/hg38/Homo_sapiens_assembly38.fasta --ref-version hg38 --data-sources-path /work/sdukoldby/resources/exome_annotation --output ECV2-4_somatic_mutect2_filtered_funcotated.vcf --output-file-format VCF
# gatk Funcotator --variant ECV2-4_somatic_mutect2_filtered.vcf --reference /work/sdukoldby/resources/hg38/Homo_sapiens_assembly38.fasta --ref-version hg38 --data-sources-path /work/sdukoldby/resources/funcotator_dataSources.v1.6.20190124s --output ECV2-4_somatic_mutect2_filtered_funcotated_full.vcf --output-file-format VCF
# gatk Funcotator --variant ECV2-4_somatic_mutect2_filtered.vcf --reference /work/sdukoldby/resources/hg38/Homo_sapiens_assembly38.fasta --ref-version hg38 --data-sources-path /work/sdukoldby/resources/funcotator_dataSources.v1.6.20190124s --output ECV2-4_somatic_mutect2_filtered_funcotated_full_gatk-update.vcf --output-file-format VCF
# gatk Funcotator --variant ECV2-4_somatic_mutect2_filtered.vcf --reference /work/sdukoldby/resources/hg38/Homo_sapiens_assembly38.fasta --ref-version hg38 --data-sources-path /work/sdukoldby/resources/funcotator_dataSources.v1.6.20190124s --output ECV2-4_somatic_mutect2_filtered_funcotated_full_dbsnp-reformat.vcf --output-file-format VCF

### dbSNP-vcf-filen for hg38 har kromosomer angivet på formatet '1' i stedet for 'chr1' som i alle andre database-filer og variant-filter
### Har derfor reformateret dbSNPfilen "hg38_All_20170710.vcf.gz" til dette format og erstattet filen i data source mappen
### Command:
### bcftools view hg38_All_20170710.vcf.gz | awk '{if($0 !~ /^#/) print "chr"$0; else print $0}' > newfile
### herefter zipped med bgzip og indexeret med tabix
