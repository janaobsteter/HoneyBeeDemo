include: multiple_genome_alignment.smk
include: qc_vcf.smk

configfile: "config/ancestral.yaml"


rule_all:
    input:
        "AncestralAllele_FinalVcf.txt"

rule rename_alleles:
    input:
        "OutspeciesInfo_All_aligned.txt"
    output:
        "AlignedSnps_focal_CarnicaChr.txt"
    script:
        "scripts"

rule sort_alleles:
    input:
        "AlignedSnps_focal_CarnicaChr.txt"
    output:
        "AlignedSnps_focal_CarnicaChrSorted.txt"
    shell:
        "sort -V {input}.txt | awk -F" " '{print $0 FS $1"_"$2}' > {output}"


rule extract_aligned_alleles_from_vcf:
    input:
        alignedAlleles="AlignedSnps_focal_CarnicaChrSorted.txt"
        vcf=config['rawVcf']
        out="RawVcfInfo"
    output:
        "RawVcfInfo.INFO"
    script:
        "scripts/ExtractPosVcf.sh {input.vcf} {input.alignedPos} {input.out}"

rule extract_pos_from_aligned_alleles_from_vcf:
    input:
        "RawVcfInfo.INFO"
    output:
        "RawVcfFullPos.txt"
    shell:
        "awk '{print $1 "_" $2}' {input} > {output}"

rule extract_vcf_alleles_from_aligned:
    input:
        alignedAlleles:"AlignedSnps_focal_CarnicaChrSorted.txt"
        vcfPos:"RawVcfFullPos.txt"
    output:
        "AlignedSnps_focal_CarnicaChrSortedRawVcf.txt"
    shell:
        "grep -Fwf {input.vcfPos} {input.alignedPop} > {output}"

rule create_estsfs_dicts:
    input:
        alignedAlleles: "AlignedSnps_focal_CarnicaChrSortedRawVcf.txt"
        vcfAlleles: "RawVcfInfo.INFO"
    output:
        "EstSfs_Dict{chunk}.csv"
    script:
        "scripts/CreateInputForEstsfs_fromWGAbed_Cactus.py"
        #CreateInputForEstsfs_Loop.sh This is qsub

rule edit_estsfs_dicts:
    input:
        "EstSfs_Dict{chunk}.csv"
    output:
        "EstSf_Dict{chunk}E.csv"
    params:
        chunk=config['noEstSfsChunks']
    shell:
        """
        for cycle in $(seq {params.chunk})
        do
          	cut -f2,3,4,5 EstSfs_Dict${cycle}.csv > EstSfs_Dict${cycle}E.csv
            grep -v "()" EstSfs_Dict${cycle}E.csv > tmp && mv tmp EstSfs_Dict${cycle}E.csv
            sed -i "s/ //g" EstSfs_Dict${cycle}E.csv
            # Set the correct separators for the file
            awk -F "\t" '{print $1"\t"$2" "$3" "$4}' EstSfs_Dict${cycle}E.csv > tmp && mv tmp EstSfs_Dict${cycle}E.csv
            # Remove the parenthesis from the file
            sed -i "s/(//g" EstSfs_Dict${cycle}E.csv
            sed -i "s/)/ /g" EstSfs_Dict${cycle}E.csv
        done
        """

rule combine_estsfs_dicts:
    input:
        expand("EstSf_Dict{chunk}E.csv", chunk = range(config['noEstSfsChunks']))
    output:
        "EstSfs_Dict_Total.csv"
    shell:
        "cat {input} > {output}"

rule run_estsfs:
    input:
        config="config-kimura_3o.txt",
        dict="EstSfs_Dict_Total.csv",
        seed="seedfile.txt"
    output:
        text="EstsfsOutput/outputEtsfs.txt",
        pvalue="EstsfsOutput/output-pvalues.txt"
    shell:
        """
        mkdir EstsfsOutput
        ./estsfs {input.config} {input.dict} {input.seed} {output.text} {output.pvalue}
        """

rule determine_ancestral_allele:
    input:
        "EstsfsOutput/output-pvalues.txt"
    output:
        "EstsfsOutput/AncestralAllele_3o.csv"
    script:
        "scripts/DetermineAncestralAllele.py"

rule edit_ancestral_allele:
    input:
        alleles:"EstsfsOutput/AncestralAllele_3o.csv"
    output:
        "EstsfsOutput/AncestralAlleles1_16.csv"
    shell:
        """
        awk -F"," '$2 != ""' {input.alleles} > tmp
        sed "s/carnica.LG//g" tmp > {output}
        """"

rule edit_ancestral_allele_positions:
    input:
        pos:"AncAllelesPos.csv"
    output:
        "EstsfsOutput/AncAllelesPos1_16.csv"
    shell:
        """
        sed "s/carnica.LG//g" {input.pos} > tmp
        awk '{print $1"_"$2}' tmp > {output}
        """

rule get_final_vcf_pos:
    input:
        vcf={config.finalVcf}
    output:
        "FinalVcf_Pos.txt"
    shell:
        """
        bcftools query -f '%CHROM %POS\n' {input}" > tmp
        awk '{print $1"_"$2}' tmp > {output}
        """

rule extract_vcf_ancAl:
    input:
        vcfPos="FinalVcf_Pos.txt"
        ancAlleles="EstsfsOutput/AncestralAlleles1_16.csv"
    output:
        "AncestralAllele_FinalVcf.txt"
    shell:
        "grep -Fwf {input.vcfPos} {input.ancAlleles} > {output}"
