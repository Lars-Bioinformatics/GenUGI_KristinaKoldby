import gzip, re

# PATH,SAMPLES,TOOLS, = glob_wildcards("{path}/{sample}.{tool}.vcf.gz")

SAMPLES = "G45-ECV2-29"
print(SAMPLES)
print(len(SAMPLES))

# TOOLS = ["manta", "lumpy"]
TOOLS = ["lumpy", "manta"]
# TOOLS_LIST = list(set(TOOLS))
# print(TOOLS_LIST)

rule all:
    input:
        #expand("survivor/{sample}.{tool}_prepared.vcf", sample=SAMPLES, tool=TOOLS),
        # expand("survivor/{sample}_merged.vcf", sample=SAMPLES)
        expand("survivor/{sample}.delly_prepared.vcf", sample=SAMPLES)

rule prepare_manta_vcf:
    input:
        vcf="manta/{sample}.manta.vcf.gz"
    output:
        vcf="survivor/{sample}.manta_prepared.vcf"
    run:
        with open(output.vcf, "w") as o:
            with gzip.open(input.vcf, "rt") as f:
                for line in f:
                    # Keep all header lines in vcf
                    if line.startswith("##"):
                        o.write(line)
                        continue
                    # Keep only tumor sample
                    lst = line.strip().split("\t")
                    if line.startswith("#CHROM"):
                        o.write("\t".join(lst[0:9]))
                        o.write("\t"+wildcards.sample+"_manta\n")
                        continue
                    # Keep only the smaller entry for translocations
                    # e.g. chr1->chr2, and not chr2->chr1
                    # ID_field = lst[2].split(":")
                    # if ID_field[0] != "MantaBND" or ID_field[7] != "0":
                    o.write("\t".join(lst[0:9]))
                    o.write("\t"+lst[10]+"\n")


rule prepare_lumpy_vcf:
    input:
        vcf="lumpy/{sample}.lumpy.vcf.gz"
    output:
        vcf="survivor/{sample}.lumpy_prepared.vcf"
    run:
        with open(output.vcf, "w") as o:
            with gzip.open(input.vcf, "rt") as f:
                for line in f:
                    # Keep all header lines in vcf
                    if line.startswith("##"):
                        o.write(line)
                        continue
                    # Keep only tumor sample
                    lst = line.strip().split("\t")
                    if line.startswith("#CHROM"):
                        o.write("\t".join(lst[0:9]))
                        o.write("\t"+wildcards.sample+"_lumpy\n")
                        continue
                    # Keep only the smaller entry for translocations
                    # e.g. chr1->chr2, and not chr2->chr1
                    # ID_field = lst[2]
                    # if ID_field.find("_2")<0:
                    o.write("\t".join(lst[0:9]))
                    o.write("\t"+lst[9]+"\n")


rule prepare_gridss_vcf:
    input:
        vcf="gridss/{sample}.gridss.vcf.gz"
    output:
        vcf="survivor/{sample}.gridss_prepared.vcf"
    run:
        with open(output.vcf, "w") as o:
            with gzip.open(input.vcf, "rt") as f:
                for line in f:
                    # Keep all header lines in vcf
                    if line.startswith("##"):
                        o.write(line)
                        continue
                    # Keep only tumor sample
                    lst = line.strip().split("\t")
                    if line.startswith("#CHROM"):
                        o.write("\t".join(lst[0:9]))
                        o.write("\t"+wildcards.sample+"_gridss\n")
                        continue
                    # Keep only the smaller entry for translocations
                    # e.g. chr1->chr2, and not chr2->chr1
                    # ID_field = lst[2]
                    # if ID_field.find("_2")<0:
                    o.write("\t".join(lst[0:9]))
                    o.write("\t"+lst[10]+"\n")

rule prepare_brass_vcf:
    input:
        vcf="brass/{sample}.brass.vcf.gz"
    output:
        vcf="survivor/{sample}.brass_prepared.vcf"
    run:
        with open(output.vcf, "w") as o:
            with gzip.open(input.vcf, "rt") as f:
                for line in f:
                    # Keep all header lines in vcf
                    if line.startswith("##"):
                        o.write(line)
                        continue
                    # Keep only tumor sample
                    lst = line.strip().split("\t")
                    if line.startswith("#CHROM"):
                        o.write("\t".join(lst[0:9]))
                        o.write("\t"+wildcards.sample+"_brass\n")
                        continue
                    # Keep only the smaller entry for translocations
                    # e.g. chr1->chr2, and not chr2->chr1
                    # ID_field = lst[2]
                    # if ID_field.find("_2")<0:
                    o.write("\t".join(lst[0:9]))
                    o.write("\t"+lst[10]+"\n")


rule prepare_delly_vcf:
    input:
        vcf="delly/{sample}.delly.vcf.gz"
    output:
        vcf="survivor/{sample}.delly_prepared.vcf"
    run:
        with open(output.vcf, "w") as o:
            with gzip.open(input.vcf, "rt") as f:
                for line in f:
                    # Keep all header lines in vcf
                    if line.startswith("##"):
                        o.write(line)
                        continue
                    # Keep only tumor sample
                    lst = line.strip().split("\t")
                    if line.startswith("#CHROM"):
                        o.write("\t".join(lst[0:9]))
                        o.write("\t"+wildcards.sample+"_delly\n")
                        continue
                    # Keep only the smaller entry for translocations
                    # e.g. chr1->chr2, and not chr2->chr1
                    if lst[2].startswith("BND"):
                        chr1 = lst[0]
                        pos1 = lst[1]
                        ref1 = lst[3]
                        
                        ## Write first break end
                        # print(lst[4].find("["))
                        if lst[4].find("[")>=0:
                            delim = "["
                        else: # ]
                            delim = "]"

                        # print(delim)
                        alt = re.split("\\"+delim+"|:", lst[4])
                        # if REF = [chrA:pos[N or ]chrA:pos]N
                        if lst[4].startswith(delim):
                            # print(alt[1:])
                            chr2, pos2, ref2 = alt[1:]
                            o.write(chr2+"\t"+pos2+"\t"+lst[2]+"_1\t"+ref2+"\t")
                            o.write(delim+chr1+":"+pos1+delim+ref1+"\t")
                        # if REF = N[chrA:pos[ or N]chrA:pos]
                        else:
                            # print(alt[:3])
                            ref2, chr2, pos2 = alt[:3]
                            o.write(chr2+"\t"+pos2+"\t"+lst[2]+"_1\t"+ref2+"\t")
                            o.write(ref1+delim+chr1+":"+pos1+delim+"\t")

                        o.write("\t".join(lst[5:9]))
                        o.write("\t"+lst[9]+"\n")
                        
                        # Write second break end - renaming ID
                        o.write("\t".join(lst[0:2]))
                        o.write("\t"+lst[2]+"_2\t") # ID
                        o.write("\t".join(lst[3:9]))
                        o.write("\t"+lst[9]+"\n")
                        
                    else:
                        # Write entire line, only selecting tumor sample
                        o.write("\t".join(lst[0:9]))
                        o.write("\t"+lst[9]+"\n")


rule merge_sv_vcf:
    input:
        vcf=expand("survivor/{{sample}}.{tool}_prepared.vcf", tool=TOOLS)
    output:
        vcf="survivor/{sample}_merged.vcf"
    shell:
        """
        ls survivor/{wildcards.sample}*vcf > survivor/{wildcards.sample}_filenames.txt
        
        SURVIVOR merge survivor/{wildcards.sample}_filenames.txt 200 0 0 0 0 0 {output.vcf}
        """

