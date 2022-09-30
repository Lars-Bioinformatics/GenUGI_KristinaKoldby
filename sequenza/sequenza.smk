
# SAMPLES, = glob_wildcards("{sample}_tumor_tagseq-medexome.connor.recalibrated.bam")
SAMPLES, = glob_wildcards("{sample}_nimblegen-medexome_HGVM5DSXX.connor.recalibrated.bam")
# SAMPLES = "ECV2-4-biopsi-H2"


configfile: "somatic_matched_samples.yaml"
ref = "/work/sdukoldby/resources/hg38/Homo_sapiens_assembly38.fasta"
output_sequenza = "sequenza/"
gcwindow = 50
seqz_bin_window = 50

# import pathlib
#
# indir = pathlib.Path("input")
# paths = indir.glob("*-?-?_?_tagseq-medexome.connor.recalibrated.bam")
# patients = set([x.stem.split("-")[1] for x in paths])
#
# def find_normal(wildcards):
#     normals = [str(x) for x in indir.glob(f"{wildcards.patients}-blod_normal_tagseq-medexome.connor.recalibrated.bam")]
#     return sorted(normals)

rule all:
    input:
        output_sequenza+"Homo_sapiens_assembly38.gc50Base.wig.gz",
        expand(output_sequenza+"{sample}.bins.seqz.gz", sample = SAMPLES),
        expand(output_sequenza+"{sample}/{sample}_gc_plots.pdf", sample = SAMPLES)

rule gc_wiggle:
    input:
        fasta = {ref}
    output:
        output_sequenza+"Homo_sapiens_assembly38.gc50Base.wig.gz"
    params:
        window = {gcwindow}
    shell:
        """
        sequenza-utils gc_wiggle \
        -w {params.window} \
        -f {input} \
        -o {output}
        """


### Process BAM and Wiggle files to produce a seqz file ###
# def get_normal(s):
#     print(s)
#     return("-".join(s.split("-")[:2])+"-blod_normal_tagseq-medexome.connor.recalibrated.bam")

rule bam2seqz:
    input:
        normalbam = "G88-F9B1-Blod_nimblegen-medexome_HGVM5DSXX.connor.recalibrated.bam",
        tumorbam = "{sample}_nimblegen-medexome_HGVM5DSXX.connor.recalibrated.bam",
        wig=output_sequenza+"Homo_sapiens_assembly38.gc50Base.wig.gz"
    output:
        seq = output_sequenza+"{sample}.seqz.gz"
    shell:
        """
        sequenza-utils bam2seqz \
        -n {input.normalbam} \
        -t {input.tumorbam} \
        --fasta {ref} \
        -gc {input.wig} \
        -o {output}
        """

### G45 exomes:
# normalbam = lambda wildcards: "-".join(wildcards.sample.split("-")[:2])+"-blod_normal_tagseq-medexome.connor.recalibrated.bam",
# tumorbam = "{sample}_tumor_tagseq-medexome.connor.recalibrated.bam",

### Post-process by binning the original seqz file ###
rule seqz_binning:
    input:
        seqz = output_sequenza+"{sample}.seqz.gz"
    output:
        bins = output_sequenza+"{sample}.bins.seqz.gz"
    params:
        window = {seqz_bin_window}
    shell:
        """
        sequenza-utils seqz_binning \
        --seqz {input} \
        -w {params.window} \
        -o {output}
        """

rule R_commands:
    input:
        output_sequenza+"{sample}.bins.seqz.gz"
    output:
        output_sequenza+"{sample}/{sample}_gc_plots.pdf"
    shell:
        """
        Rscript --vanilla sequenza/sequenza.R \
        {input} \
        {wildcards.sample}
        """
