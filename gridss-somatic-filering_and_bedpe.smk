from snakemake.utils import R

PON = "small_pon"

# SAMPLES, = glob_wildcards("{sample}__"+PON+".somatic.annotated.vcf")
SAMPLES, = glob_wildcards("{sample}.vcf.gz")
print(SAMPLES)

REF = "/work/sdukoldby/resources/hg38/hg38_new_names/Homo_sapiens_assembly38.fasta"
SCRIPTS = "/work/sdukoldby/data/G45-2016_genugi/190325_A00653_0016_AHGL2LDSXX/BaseCalls/hg38_withoutMinusM/GRIDSS_joint_calling/gridss-2.8.1/"

mem = 12

rule all:
    input:
        # expand("{sample}__"+PON+".filtered.somatic.vcf.gz", sample=SAMPLES)
        expand("{sample}_gridss_single_breakends.bedpe", sample=SAMPLES),
        # expand("{sample}__"+PON+".filtered.somatic.annotated.vcf.gz", sample=SAMPLES)

# rule remove_germline:
#     input:
#         "{sample}__"+PON+".filtered.vcf.gz"
#     output:
#         "{sample}__"+PON+".filtered.somatic.vcf.gz"
#     shell:
#         """
#         bcftools filter \
#         --include 'FORMAT/VF[0:0] = 0' \
#         {input} | \
#         bgzip -c > {output}
#         """

# rule annotateUntemplated_variants:
#     input:
#         vcf="{sample}__"+PON+".filtered.somatic.vcf.gz"
#     output:
#         vcf="{sample}__"+PON+".filtered.somatic.annotated.vcf.gz"
#     threads: 24
#     shell:
#         """
#          java -Xmx{mem}g \
#             -cp {SCRIPTS}/gridss-2.8.1-gridss-jar-with-dependencies.jar \
#             gridss.AnnotateUntemplatedSequence \
#             REFERENCE_SEQUENCE={REF} \
#             INPUT={input} \
#             OUTPUT={output} \
#             WORKER_THREADS={threads}
#         """

# rule prepare_reference_for_repeatMasker_annotation:

# rule repeatMasker_annotation:
#     input:
#         gridss_ref="Homo_sapiens_assembly38.fasta.out",
#         vcf="{sample}__"+PON+".filtered.somatic.vcf.gz"
#     output:
#         vcf="{sample}__"+PON+".filtered.somatic.annotated.vcf.gz"
#     shell:
#         """
#         Rscript {SCRIPTS}/gridss_annotate_insertions_repeatmaster.R \
#             --input {input.vcf} \
#             --output {output.vcf} \
#             --repeatmasker {input.gridss_ref} \
#             --scriptdir {SCRIPTS}
#         """

rule create_bedpe_and_bed:
    input:
        # vcf="{sample}__"+PON+".somatic.annotated.vcf"
        vcf="{sample}.vcf.gz"
    output:
        bedpe="{sample}_gridss_breakpoints.bedpe",
        bed="{sample}_gridss_single_breakends.bed"
    run:
        # From: https://github.com/PapenfussLab/gridss/blob/master/example/bedpe-export-example.R
        R("""
        library(StructuralVariantAnnotation)
        library(rtracklayer)

        vcf = readVcf("{input}")

        # Export breakpoints to BEDPE
        bpgr = breakpointRanges(vcf)
        # TODO: add your event filtering here. The default GRIDSS output is very verbose/sensitive.
        write.table(breakpointgr2bedpe(bpgr), file="{output.bedpe}", sep="\\t", quote=FALSE, col.names=FALSE, row.names=FALSE)

        # Export single breakends to BED
        begr = breakendRanges(vcf)
        # TODO: add your event filtering here. The default GRIDSS output is very verbose/sensitive.
        begr$score = begr$QUAL
        export(begr, con="{output.bed}")
        """)
    # script:
    #     "scripts/create_bedpe_and_bed_files.R"

rule single_breaks_to_bedpe:
    input:
        bed="{sample}_gridss_single_breakends.bed"
    output:
        bedpe="{sample}_gridss_single_breakends.bedpe"
    shell:
        """
        awk -F'\\t' '{{ \
            printf "%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n", $1, $2, $2+1, $1, $3-1, $3, $4, $5, $6, $6 \
        }}' {input} > {output}
        """
