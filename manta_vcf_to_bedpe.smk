SAMPLES, = glob_wildcards("{sample}.vcf.gz")
# SAMPLES = "G45-ECV2-31"

rule all:
    input:
        expand("bedpe/{sample}.bedpe", sample=SAMPLES),
        expand("bedpe/{sample}.deletion.bedpe", sample=SAMPLES)
        
rule convertToBedpe:
    input:
        vcf="{sample}.vcf.gz"
    output:
        bedpe="bedpe/{sample}.bedpe"
    shell:
        """
        mkdir -p bedpe
        
        zcat {input.vcf} > {wildcards.sample}.vcf
        SURVIVOR vcftobed {wildcards.sample}.vcf 1 10000000000000000 {wildcards.sample}.bed
        
        # Fix bed coordinates - and keep only one entry for translocations
        cat {wildcards.sample}.bed | \
        awk '{{split($7,a,":"); if (!($11=="TRA" && a[8]==1)) print $0}}' | \
        awk -F'\\t' '{{for (i=1;i<=NF;i++){{
            if (i==2 || i==5)       printf "%s\\t", $i-2;
            else if (i==3 || i==6)  printf "%s\\t", $i-1;
            else if (i==NF)         printf "%s\\n", $i;
            else                    printf "%s\\t", $i;
        }} }}' > {output.bedpe}
    
        rm -f {wildcards.sample}.bed {wildcards.sample}.vcf
        """

rule seperate_sv_type:
    input:
        bedpe="bedpe/{sample}.bedpe"
    output:
        deletions="bedpe/{sample}.deletion.bedpe",
        inversions="bedpe/{sample}.inversion.bedpe",
        translocations="bedpe/{sample}.translocation.bedpe",
        duplications="bedpe/{sample}.duplication.bedpe",
        insertions="bedpe/{sample}.insertion.bedpe"
    shell:
        """
        cat {input} | grep 'DEL' | cut -f 1-6 | sed 's/chr/hs/g' > {output.deletions};
        cat {input} | grep 'INV' | cut -f 1-6 | sed 's/chr/hs/g' > {output.inversions};
        cat {input} | grep 'TRA' | cut -f 1-6 | sed 's/chr/hs/g' > {output.translocations};
        cat {input} | grep 'DUP' | cut -f 1-6 | sed 's/chr/hs/g' > {output.duplications};
        cat {input} | grep 'INS' | cut -f 1-6 | sed 's/chr/hs/g' > {output.insertions}
        """