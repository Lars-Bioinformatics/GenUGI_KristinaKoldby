import gzip, re, sys
from snakemake.utils import R
from collections import Counter

PATH,SAMPLES,TOOLS, = glob_wildcards("{path}/{sample}.{tool}.vcf.gz")
SAMPLES = sorted(list(set(SAMPLES)))

# SAMPLES = "G45-ECV2-4"
# SAMPLES = "G45-ECV2-31"
print(SAMPLES)
print(len(SAMPLES))

# TOOLS = ["delly", "gridss", "lumpy", "manta"]
# TOOLS = ["brass", "lumpy", "manta"]
# TOOLS = ["lumpy", "manta"]
TOOLS = sorted(list(set(TOOLS)))
print(TOOLS)

# Tool pairs
PAIRS = [(TOOLS[i],TOOLS[j]) for i in range(len(TOOLS)) for j in range(i+1,len(TOOLS))]
print(PAIRS)

# Slop (i.e. distance) around breakends to include in overlap
SLOP="200"

onstart:
    shell("mkdir -p survivor survivor/{{prepared_bedpe,pairwise_merged_bedpe,unique_SVs}}")


###############################################################################
### Rules
###############################################################################
rule all:
    input:
        expand("survivor/{sample}.merged.SLOP_{slop}.txt", sample=SAMPLES, slop=SLOP),
        expand("survivor/unique_SVs/{sample}_unique_brass.SLOP_{slop}.txt", sample=SAMPLES, slop=SLOP),
        expand("survivor/unique_SVs/{sample}_unique_manta.SLOP_{slop}.txt", sample=SAMPLES, slop=SLOP),
        expand("survivor/unique_SVs/{sample}_only_in_{tool}.SLOP_{slop}.txt", sample=SAMPLES, slop=SLOP, tool=["delly","lumpy","gridss"])
        # expand("survivor/prepared_bedpe/{sample}.{tool}_prepared.bedpe", sample=SAMPLES, tool=TOOLS),
        # [expand("survivor/pairwise_merged_bedpe/{sample}.{tool1}_and_{tool2}.SLOP_"+SLOP+".bedpe",
        #     sample=SAMPLES, tool1=t1, tool2=t2) for (t1,t2) in PAIRS]
        

###############################################################################
### Convert vcf to bedpe
###############################################################################
rule prepare_manta_bedpe:
    input:
        vcf="manta/{sample}.manta.vcf.gz"
    output:
        bedpe="survivor/prepared_bedpe/{sample}.manta_prepared.bedpe"
    shell:
        """
        # Filter PASS variants and convert vcf to bedpe
        zcat {input.vcf} | grep PASS > {wildcards.sample}_manta_tmp.vcf
        SURVIVOR vcftobed {wildcards.sample}_manta_tmp.vcf 1 10000000000000000 {wildcards.sample}.bed
        
        # Fix bed coordinates - and keep only one entry for translocations
        awk '{{split($7,a,":"); if (!($11=="TRA" && a[8]==1)) print $0}}' {wildcards.sample}.bed | \
        awk -F'\\t' '{{for (i=1;i<=NF;i++){{
            if (i==2 || i==5)       printf "%s\\t", $i-2;
            else if (i==3 || i==6)  printf "%s\\t", $i-1;
            else if (i==NF)         printf "%s\\n", $i;
            else                    printf "%s\\t", $i;
        }} }}' > {output.bedpe}
        
        # remove tmp files
        rm -f {wildcards.sample}.bed {wildcards.sample}_manta_tmp.vcf
        """


rule prepare_lumpy_bedpe:
    input:
        vcf="lumpy/{sample}.lumpy.vcf.gz"
    output:
        bedpe="survivor/prepared_bedpe/{sample}.lumpy_prepared.bedpe"
    shell:
        """
        # Convert vcf to bedpe
        zcat {input.vcf} > {wildcards.sample}_lumpy_tmp.vcf
        SURVIVOR vcftobed {wildcards.sample}_lumpy_tmp.vcf 1 10000000000000000 {wildcards.sample}.bed
        
        # Fix bed coordinates - and keep only one entry for translocations
        grep -v _2 {wildcards.sample}.bed | \
        awk -F'\\t' '{{for (i=1;i<=NF;i++){{
            if (i==2 || i==5)       printf "%s\\t", $i-2;
            else if (i==3 || i==6)  printf "%s\\t", $i-1;
            else if (i==NF)         printf "%s\\n", $i;
            else                    printf "%s\\t", $i;
        }} }}' > {output.bedpe}
    
        # remove tmp files
        rm -f {wildcards.sample}.bed {wildcards.sample}_lumpy_tmp.vcf
        """
    # else if (i==7)          printf "%s\\t", gsub("_1","",$i);
    # else if (i==8)         continue; # skip commas


rule prepare_gridss_bedpe:
    input:
        vcf="gridss/{sample}.gridss.vcf.gz"
    output:
        bedpe="survivor/prepared_bedpe/{sample}.gridss_prepared.bedpe"
    run:
        R("""
        library(StructuralVariantAnnotation)
        library(rtracklayer)

        vcf = readVcf("{input.vcf}")

        # Export breakpoints to BEDPE
        bpgr = breakpointRanges(vcf)
        bedpe = breakpointgr2bedpe(bpgr)
        bedpe$svclass = info(vcf)[bedpe$name,]$SIMPLE_TYPE

        # Export single breakends to BED
        begr2 = breakendRanges(vcf)
        begr2$score = begr2$QUAL
        export(begr2, con="{wildcards.sample}_gridss_tmp.bed")
        bed=read.table("{wildcards.sample}_gridss_tmp.bed", sep = "\\t")
        bedpe_unk=data.frame(bed$V1,bed$V2,bed$V3,bed$V1,bed$V2,bed$V3+begr2$insLen,bed$V4,bed$V5,bed$V6,bed$V6,"UNK")
        names(bedpe_unk) = names(bedpe)

        # Merge SVs
        bedpe_final = rbind(bedpe, bedpe_unk)
        bedpe_final$svclass[bedpe_final$svclass=="CTX"] = "TRA" # rename translocations
        write.table(bedpe_final, file="{output.bedpe}", sep="\\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
        
        # Remove tmp file
        file.remove("{wildcards.sample}_gridss_tmp.bed")
        """)

rule prepare_brass_bedpe:
    input:
        bedpe="brass/bedpe/{sample}.brass.annot.bedpe.gz"
    output:
        bedpe="survivor/prepared_bedpe/{sample}.brass_prepared.bedpe"
    shell:
        """
        # Select needed columns, sort chr numbers numerically and rename SV class
        zcat {input.bedpe} | grep -v ^# | cut -f 1-10,12 | sort -k1,1 -k2,2n -V | \
        sed -e 's/translocation/TRA/g' -e 's/deletion/DEL/g' \
            -e 's/inversion/INV/g' -e 's/tandem-duplication/DUP/g'> {output.bedpe}
        """


rule prepare_delly_bedpe:
    input:
        vcf="delly/{sample}.delly.vcf.gz"
    output:
        bedpe="survivor/prepared_bedpe/{sample}.delly_prepared.bedpe"
    shell:
        """
        zcat {input.vcf} > {wildcards.sample}_delly_tmp.vcf
        SURVIVOR vcftobed {wildcards.sample}_delly_tmp.vcf 1 10000000000000000 {wildcards.sample}.bed
        
        # Fix bed coordinates - and keep only one entry for translocations
        cat {wildcards.sample}.bed | \
        awk -F'\\t' '{{for (i=1;i<=NF;i++){{
            if (i==2 || i==5)       printf "%s\\t", $i-2;
            else if (i==3 || i==6)  printf "%s\\t", $i-1;
            else if (i==NF)         printf "%s\\n", $i;
            else                    printf "%s\\t", $i;
        }} }}' > {output.bedpe}
        
        # remove tmp files
        rm -f {wildcards.sample}.bed {wildcards.sample}_delly_tmp.vcf
        """

###############################################################################
### Merge SV calls
###############################################################################
rule pairwise_merge_bedpe_SVs:
    input:
        bedpe1="survivor/prepared_bedpe/{sample}.{tool1}_prepared.bedpe",
        bedpe2="survivor/prepared_bedpe/{sample}.{tool2}_prepared.bedpe"
    output:
        bedpe="survivor/pairwise_merged_bedpe/{sample}.{tool1}_and_{tool2}.SLOP_{slop}.bedpe"
    shell:
        # bedtools pairtopair (aka pairToPair)
        # -is: no strand
        # -slop N: allowed bases between breakends
        # -type both: SVs overlaps both breakends
        # $11==$22: check SV class identical
        """
        pairToPair -slop {SLOP} -is -type both -a {input.bedpe1} -b {input.bedpe2} | \
        awk '{{if ($11==$22) print}}' > {output.bedpe}
        """

# rule pairwise_merge_bedpe_SVs:
#     input:
#         bedpe1="survivor/prepared_bedpe/{sample}.{tool1}_prepared.bedpe",
#         bedpe2="survivor/prepared_bedpe/{sample}.{tool2}_prepared.bedpe"
#     output:
#         bedpe="survivor/pairwise_merged_bedpe/{sample}.{tool1}_and_{tool2}.SLOP_{slop}.bedpe"
#     shell:
#         # bedtools pairtopair (aka pairToPair)
#         # -is: no strand
#         # -slop N: allowed bases between breakends
#         # -type both: SVs overlaps both breakends
#         # $11==$22: check SV class identical
#         """
#         pairToPair -slop {SLOP} -is -type either -a {input.bedpe1} -b {input.bedpe2} > {output.bedpe}
#         """


rule merge_all_bedpe_SVs:
    input:
        [expand("survivor/pairwise_merged_bedpe/{{sample}}.{tool1}_and_{tool2}.SLOP_{{slop}}.bedpe",
            tool1=t1, tool2=t2) for (t1,t2) in PAIRS]
    output:
        txt="survivor/{sample}.merged.SLOP_{slop}.txt"
    run:
        def traverse_graph(graph, cur_node, visited):
            # Add current node to list of visited nodes
            visited.add(cur_node)
            # Visit each node connected to current node 
            for node in graph[cur_node]:
                if node not in visited:
                    graph, visited = traverse_graph(graph, node, visited)
                    del graph[node]
            return (graph, visited)


        # Build graph from pairwise merged SV calls
        graph = dict()
        for i,in_file in enumerate(input):
            # print(in_file)
            with open(in_file, "r") as f:
                t1, t2 = PAIRS[i]
                # print(t1,t2)
                for line in f:
                    lst = line.strip().split("\t")
                    # print(lst)
                    # SV_id, caller_tool, sv_class
                    id1 = (lst[6],t1,lst[10])
                    id2 = (lst[17],t2,lst[21])
                    # Add id2 to id1's list of links
                    if id1 in graph:
                        graph[id1].add(id2)
                    else:
                        graph[id1] = set()
                        graph[id1].add(id2)
                    
                    # Add id1 to id2's list of links
                    if id2 in graph:
                        graph[id2].add(id1)
                    else:
                        graph[id2] = set()
                        graph[id2].add(id1)


        nodes = [key for key in graph]
        print(len(nodes))
        # print(graph)
        o = open(output.txt, "w")
        o.write("Type\t" + '\t'.join(TOOLS) + "\n")
        for node in nodes:
            if node in graph:
                graph, visited = traverse_graph(graph, node, set())
                if len(visited)>0:
                    # Create info object to write tool ids to file
                    info = {"type":list(visited)[0][2]}
                    for v in visited:
                        # v[1] = tool 
                        if v[1] in info:
                            info[v[1]].append(v[0])
                        else:
                            info[v[1]] = [v[0]]

                    tool_ids = [','.join(info.get(tool, ["NA"])) for tool in TOOLS]

                    print([info["type"]] + tool_ids)
                    o.write('\t'.join([info["type"]] + tool_ids))
                    o.write('\n')   
        o.close()
            

###############################################################################
### Get list of SVs unique to each tool
###############################################################################
rule get_unique_breaks_brass:
    input:
        txt=expand("survivor/{{sample}}.merged.SLOP_{slop}.txt", slop=SLOP),
        bedpe="survivor/prepared_bedpe/{sample}.brass_prepared.bedpe",
        vcf="brass/{sample}.brass.vcf.gz"
    output:
        txt="survivor/unique_SVs/{sample}_unique_brass.SLOP_{slop}.txt"
    run:
        # Get list of IDs also found in other tools
        with open(input[0], "r") as f:
            header = f.readline()
            tool_index = header.split("\t").index("brass")
            # print(tool_index)
            nonUniqueSVs = set()
            for line in f:
                # print(line.split("\t")[tool_index].split(","))
                for e in line.split("\t")[tool_index].split(","):
                    if e != "NA":
                        nonUniqueSVs.add(e)
            # print(nonUniqueSVs)
        
        # Get list of IDs having a balanced translocation (we only want to count once)
        def find_between( s, first, last ):
            try:
                start = s.index( first ) + len( first )
                end = s.index( last, start )
                return s[start:end]
            except ValueError:
                return ""


        def extract_ids(s):
            val = s.split("|")
            ids = set()
            for v in val:
                if v.index("_") == v.rindex("_"):
                    ids.add(v.split("_")[0])
                else:
                    # print(v)
                    split = v.index("_") + 2
                    # print(v[:split])
                    # print(v[split:])
                    ids.add(v[:split].split("_")[0])
                    ids.add(v[split:].split("_")[0])
            return ids


        def collapse_ids(bals_ids, cur_var, collapsedSVs):
            # Visit each balanced translocation  
            for var in bals_ids[cur_var]:
                if var not in collapsedSVs:
                    # Add current variant to list of unique SVs
                    collapsedSVs.add(var)
                    if var in bals_ids:
                        bals_ids, collapsedSVs = collapse_ids(bals_ids, var, collapsedSVs)
                        del bals_ids[var]
            return (bals_ids, collapsedSVs)

        
        def isUniqueSV(SVs, nonUniqueSVs):
            for sv in SVs:
                if sv in nonUniqueSVs:
                    return(False)
            return(True)


        with open(output.txt, "w") as o:
            o.write("Type\tbrass_id\n")
            # uniqueSVs = set()

            # Find balanced translocations
            with gzip.open(input[2], "rt") as f:
                bals_ids = dict()
                for line in f:
                    if not line.startswith("#"):
                        id, index = line.split("\t")[2].split("_")
                        if index != 2:
                            res = find_between(line, "BALS=", ";")
                            if res != "":
                                ids = extract_ids(res)
                                bals_ids[id] = ids
                # print(bals_ids)
                
                # Collapse balanced translocations
                ids = [key for key in bals_ids]
                for var in ids:
                    if var in bals_ids:
                        bals_ids, collapsedSVs = collapse_ids(bals_ids, var, set([var]))
                        # print(collapsedSVs)
                        # Check that none of the ids have been found in merged set of SVs
                        if isUniqueSV(collapsedSVs, nonUniqueSVs):
                            o.write("TRA\t" + ','.join(sorted(list(collapsedSVs))) + "\n")
                            # uniqueSVs = uniqueSVs.union(collapsedSVs)
                        nonUniqueSVs.union(collapsedSVs)
                        
                # print(uniqueSVs)
                
            
            # Get list of IDs not found in other tools (except balenced translocations)
            with open(input[1], "r") as f:
                for line in f:
                    val = line.strip().split("\t")
                    if val[6] not in nonUniqueSVs:
                        # nonUniqueSVs.add(val[6])
                        o.write(val[10] + "\t" + val[6] + "\n")
                    

rule get_unique_breaks_manta:
    input:
        txt=expand("survivor/{{sample}}.merged.SLOP_{slop}.txt", slop=SLOP),
        bedpe="survivor/prepared_bedpe/{sample}.manta_prepared.bedpe"
    output:
        txt="survivor/unique_SVs/{sample}_unique_manta.SLOP_{slop}.txt"
    run:
        # Get list of IDs also found in other tools
        with open(input[0], "r") as f:
            header = f.readline()
            tool_index = header.strip().split("\t").index("manta")
            nonUniqueSVs = set()
            for line in f:
                for e in line.strip().split("\t")[tool_index].split(","):
                    if e != "NA":
                        e = ':'.join(e.split(":")[0:4])
                        nonUniqueSVs.add(e)
        # print(nonUniqueSVs)
        
        uniqueSVs = dict()
        with open(input[1], "r") as f:
            for line in f:
                val = line.strip().split("\t")
                id = ':'.join(val[6].split(":")[0:4])
                if id not in nonUniqueSVs:
                    tup = (id, val[10])
                    if tup in uniqueSVs:
                        uniqueSVs[tup].append(val[6])
                    else:
                        uniqueSVs[tup] = [val[6]]

        with open(output.txt, "w") as o:
            o.write("Type\tmanta_id\n")
            # Get list of IDs not found in other tools
            for tup in uniqueSVs:
                o.write(tup[1] + "\t" + ','.join(sorted(uniqueSVs[tup])) + "\n")


# GRIDSS, Lumpy, Delly
rule get_unique_breaks:
    input:
        txt=expand("survivor/{{sample}}.merged.SLOP_{slop}.txt", slop=SLOP),
        bedpe="survivor/prepared_bedpe/{sample}.{tool}_prepared.bedpe"
    output:
        txt="survivor/unique_SVs/{sample}_only_in_{tool}.SLOP_{slop}.txt"
    run:
        # Get list of IDs also found in other tools
        with open(input[0], "r") as f:
            header = f.readline()
            tool_index = header.split("\t").index(wildcards.tool)
            nonUniqueSVs = set()
            for line in f:
                for e in line.split("\t")[tool_index].split(","):
                    if e != "NA":
                        nonUniqueSVs.add(e)

        with open(output.txt, "w") as o:
            o.write("Type\t" + wildcards.tool + "_id\n")
            # Get list of IDs not found in other tools
            with open(input[1], "r") as f:
                for line in f:
                    val = line.strip().split("\t")
                    if val[6] not in nonUniqueSVs:
                        nonUniqueSVs.add(val[6])
                        o.write(val[10] + "\t" + val[6] + "\n")



