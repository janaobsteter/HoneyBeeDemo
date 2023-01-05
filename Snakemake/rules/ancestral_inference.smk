include: multiple_genome_alignment.smk
include: qc_vcf.smk

configfile: "config/ancestral.yaml"


rule_all:
    input:
        "AncestralAllele_FinalVcf.txt"

rule rename_alleles:
    input:
        "MultipleGenomeAlignment/OutspeciesInfo_All_aligned.txt"
    output:
        "MultipleGenomeAlignment/AlignedSnps_focal_CarnicaChr.txt"
    script:
        "scripts"

rule sort_alleles:
    input:
        "MultipleGenomeAlignment/AlignedSnps_focal_CarnicaChr.txt"
    output:
        "MultipleGenomeAlignment/AlignedSnps_focal_CarnicaChrSorted.txt"
    shell:
        "sort -V {input}.txt | awk -F" " '{print $0 FS $1"_"$2}' > {output}"


rule extract_aligned_alleles_from_vcf:
    input:
        alignedAlleles="MultipleGenomeAlignment/AlignedSnps_focal_CarnicaChrSorted.txt",
        vcf=config['rawVcf'],
        out="RawVcfInfo"
    output:
        "AncestralAllele/RawVcfInfo.INFO"
    script:
        "scripts/ExtractPosVcf.sh {input.vcf} {input.alignedPos} {input.out}"

rule extract_pos_aligned_alleles_from_vcf:
    input:
        "AncestralAllele/RawVcfInfo.INFO"
    output:
        "AncestralAllele/RawVcfFullPos.txt"
    shell:
        "awk '{print $1 "_" $2}' {input} > {output}"

rule extract_vcf_alleles_from_aligned:
    input:
        alignedAlleles="MultipleGenomeAlignment/AlignedSnps_focal_CarnicaChrSorted.txt",
        vcfPos="AncestralAllele/RawVcfFullPos.txt"
    output:
        "AncestralAllele/AlignedSnps_focal_CarnicaChrSortedRawVcf.txt"
    shell:
        "grep -Fwf {input.vcfPos} {input.alignedPop} > {output}"

rule create_estsfs_dicts:
    input:
        alignedAlleles="AncestralAllele/lignedSnps_focal_CarnicaChrSortedRawVcf.txt",
        vcfAlleles="AncestralAllele/RawVcfInfo.INFO"
    output:
        expand("AncestralAllele/Estsfs/EstSfs_Dict{{chunk}}.csv")
    script:
        "scripts/CreateInputForEstsfs_fromWGAbed_Cactus.py"
        #CreateInputForEstsfs_Loop.sh This is qsub

rule edit_estsfs_dicts:
    input:
        "AncestralAllele/Estsfs//EstSfs_Dict{chunk}.csv"
    output:
        "AncestralAllele/Estsfs//EstSf_Dict{chunk}E.csv"
    params:
        chunk=config['noEstSfsChunks']
    shell:
        """
      	cut -f2,3,4,5 {input} |
        grep -v "()" |
        sed -i "s/ //g" |
        # Set the correct separators for the file
        awk -F "\t" '{print $1"\t"$2" "$3" "$4}' |
        # Remove the parenthesis from the file
        sed -i "s/(//g" |
        sed -i "s/)/ /g" > {output}
        """

rule combine_estsfs_dicts:
    input:
        expand("AncestralAllele/Estsfs/EstSf_Dict{chunk}E.csv", chunk = range(config['noEstSfsChunks']))
    output:
        "AncestralAllele/Estsfs/EstSfs_Dict_Total.csv"
    shell:
        "cat {input} > {output}"

rule run_estsfs:
    input:
        config="AncestralAllele/Estsfs/config-kimura_3o.txt",
        dict="AncestralAllele/Estsfs/EstSfs_Dict_Total.csv",
        seed="AncestralAllele/Estsfs/seedfile.txt"
    output:
        text="AncestralAllele/Estsfs/outputEtsfs.txt",
        pvalue="AncestralAllele/Estsfs/output-pvalues.txt"
    shell:
        """
        ./estsfs {input.config} {input.dict} {input.seed} {output.text} {output.pvalue}
        """

rule determine_ancestral_allele:
    input:
        "AncestralAllele/Estsfs/output-pvalues.txt"
    output:
        "AncestralAllele/AncestralAllele_3o.csv"
    script:
        "scripts/DetermineAncestralAllele.py"

rule edit_ancestral_allele:
    input:
        alleles:"AncestralAllele/AncestralAllele_3o.csv"
    output:
        "AncestralAllele/AncestralAlleles1_16.csv"
    shell:
        """
        awk -F"," '$2 != ""' {input.alleles} > tmp
        sed "s/carnica.LG//g" tmp > {output}
        """"

rule edit_ancestral_allele_positions:
    input:
        pos:"AncestralAllele/AncAllelesPos.csv" #KJE DOBIŠ TO????
    output:
        "AncestralAllele/AncAllelesPos1_16.csv"
    shell:
        """
        sed "s/carnica.LG//g" {input.pos} > tmp
        awk '{print $1"_"$2}' tmp > {output}
        """

rule get_final_vcf_pos:
    input:
        vcf={config.finalVcf}
    output:
        "AncestralAllele/FinalVcf_Pos.txt"
    shell:
        """
        bcftools query -f '%CHROM %POS\n' {input}" > tmp
        awk '{print $1"_"$2}' tmp > {output}
        """

rule extract_vcf_ancAl:
    input:
        vcfPos="AncestralAllele/FinalVcf_Pos.txt",
        ancAlleles="AncestralAllele//AncestralAlleles1_16.csv"
    output:
        "AncestralAllele/AncestralAllele_FinalVcf.txt"
    shell:
        "grep -Fwf {input.vcfPos} {input.ancAlleles} > {output}"
