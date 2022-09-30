## Run from folder with recalibrated bam files ##
configfile: "somatic_matched_samples.yaml"
# configfile: "somatic_matched_samples_oneSample.yaml"
# configfile: "somatic_matched_samples_test.yaml"

PAIR, = glob_wildcards("{sample}_tumor_tagseq-medexome.connor.recalibrated.bam")
# PAIR = ["ECV2-8-biopsi-I1"]
# PAIR = ["ECV2-29-biopsi-C1_tumor"]
PAIR = [pair+"_tumor" for pair in PAIR]
print(PAIR)

onstart:
    shell("mkdir -p ngcgh")

regions = "/work/sdukoldby/resources/hg38/MedExome_hg38_capture_targets.bed"
window_size = 1000

rule all:
    input:
        [expand("ngcgh_lars/{tumor}_vs_{normal}_tagseq-medexome_ngcgh_w{window_size}.txt",
                    window_size=window_size,
        			normal=config[pair]["normal"],
        			tumor=config[pair]["tumor"]) for pair in PAIR]

rule ngcgh:
    input:
        normal = "{normal}_tagseq-medexome.connor.recalibrated.bam",
        tumor = "{tumor}_tagseq-medexome.connor.recalibrated.bam",
    output:
        "ngcgh_lars/{tumor}_vs_{normal}_tagseq-medexome_ngcgh_w{window_size}.txt"
    shell:
        """
        ngCGH {input.normal} {input.tumor} \
        -w {window_size} \
        -o {output}
        """
