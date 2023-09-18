from __future__ import division
import sys

configfile: "../config/ancestral_Eddie.yaml"

rule all:
    input:
        "AncestralAllele/AncestralAllele_Vcf.txt"

rule extract_aligned_alleles_from_vcf:
    input:
        alignedFocal=config['alignedFocal'],
        vcf=config['rawVcf']
    output:
        info="AncestralAllele/RawVcfInfo.INFO"
    params:
        prefix="AncestralAllele/RawVcfInfo"
    conda: "bcftools"
    shell:
        "scripts/ExtractPosVcf.sh {input.vcf} {input.alignedFocal} {params.prefix}"


rule extract_snps_from_info:
    input:
        "AncestralAllele/RawVcfInfo.INFO"
    output:
        "AncestralAllele/RawVcfInfo_SNPs.INFO"
    shell:
        """
        head -n1 {input} > Header
        sed -i "s/\t/ /g" Header
        ./scripts/ExtractSNPsFromInfo.awk {input} > tmp
        cat Header tmp > {output}
        rm Header tmp
        """

rule extract_pos_aligned_alleles_from_vcf:
    input:
        rules.extract_snps_from_info.output
    output:
        "AncestralAllele/RawVcfFullPos.txt"
    shell:
        """
        awk '{{print $1" "$2}}' {input} > {output}
        """

rule extract_vcf_alleles_from_aligned:
    input:
        alignedFocal=config['alignedFocal'],
        vcfPos=rules.extract_pos_aligned_alleles_from_vcf.output
    output:
        "AncestralAllele/AlignedSnps_focal_SortedRawVcf.txt"
    shell:
        """
        grep -Fwf {input.vcfPos} {input.alignedFocal} > {output}
        """

rule create_estsfs_dicts:
    input:
        alignedAlleles=rules.extract_vcf_alleles_from_aligned.output,
        vcfAlleles=rules.extract_snps_from_info.output
    conda: "HLab_tsinfer"
    output:
        expand("AncestralAllele/Estsfs/EstSfs_Dict{{chunk}}.csv")
    wildcard_constraints:
        chunk="\d+"
    params:
        outDir="AncestralAllele/Estsfs",
        noCycle=config['noEstSfsChunks']
    shell:
        "python scripts/CreateInputForEstsfs_fromWGAbed_Cactus.py {wildcards.chunk} {params.noCycle} {input.alignedAlleles} {input.vcfAlleles} {params.outDir}"
        #CreateInputForEstsfs_Loop.sh This is qsub

rule edit_estsfs_dicts:
    input:
        rules.create_estsfs_dicts.output
    output:
        "AncestralAllele/Estsfs/EstSfs_Dict{chunk}E.csv"
    params:
        chunk=config['noEstSfsChunks']
    wildcard_constraints:
        chunk="\d+"
    log:
        "logs/EditDict{chunk}.log"
    shell:
        """
      	cut -f2,3,4,5 {input} |	grep -v "()" > tmp1
        sed -i "s/ //g" tmp1 
        # Set the correct separators for the file
        awk -F "\t" '{{print $1"\t"$2" "$3" "$4}}' tmp1 > tmp2
        # Remove the parenthesis from the file
        sed -i "s/(//g" tmp2 
        sed "s/)/ /g" tmp2 > {output}
        rm tmp1 tmp2 
        """

rule combine_estsfs_dicts:
    input:
        expand("AncestralAllele/Estsfs/EstSfs_Dict{chunk}E.csv", chunk = range(config['noEstSfsChunks']))
    output:
        "AncestralAllele/Estsfs/EstSfs_Dict_Total.csv"
    wildcard_constraints:
        chunk="\d+"
    shell:
        "cat {input} > {output}"

rule run_estsfs:
    input:
        dict="AncestralAllele/Estsfs/EstSfs_Dict_Total.csv",
        config=config['estsfsConfig'],
        seed=config['estsfsSeed']
    output:
        text="AncestralAllele/Estsfs/outputEtsfs.txt",
        pvalue="AncestralAllele/Estsfs/output-pvalues.txt"
    shell:
        """
        ./scripts/est-sfs-release-2.04/est-sfs {input.config} {input.dict} {input.seed} {output.text} {output.pvalue}
        """

rule extract_major:
    input:
        "AncestralAllele/Estsfs/EstSfs_Dict_Total.csv"
    output:
        "AncestralAllele/MajorAllele.txt"
    shell:
        "./scripts/ExtractMajorAllele.awk {input} > {output}"

rule extract_minor:
    input:
        "AncestralAllele/Estsfs/EstSfs_Dict_Total.csv"
    output:
        "AncestralAllele/MinorAllele.txt"
    shell:
        "./scripts/ExtractMinorAllele.awk {input} > {output}"

rule extract_major_outgroup:
    input:
        "AncestralAllele/Estsfs/EstSfs_Dict_Total.csv"
    output:
        "AncestralAllele/MajorOutgroup.txt"
    shell:
        "./scripts/ExtractMajorOutgroup.awk {input} > {output}"

rule extract_estsfs_prob:
    input:
        "AncestralAllele/Estsfs/output-pvalues.txt"
    output:
        "AncestralAllele/AncestralProb.txt"
    shell:
        "tail -n +9 {input} | cut -f3 -d' ' > {output}"

rule determine_ancestral:
    input:
        ancProb="AncestralAllele/AncestralProb.txt",
        major="AncestralAllele/MajorAllele.txt",
        minor="AncestralAllele/MinorAllele.txt",
        majorOut="AncestralAllele/MajorOutgroup.txt"
    output:
        "AncestralAllele/AncestralAllele_Vcf.txt"
    shell:
        """
        paste {input.ancProb} {input.major} {input.minor} {input.majorOut} > AncestralAllele/ProbMajorMinorOutgroup.txt
        awk '{{ if(($1 >= 0.5)) {{print $2}} else if (($1 < 0.5)) {{print $3}} else {{print $4}} }}' AncestralAllele/ProbMajorMinorOutgroup.txt > {output}
        #rm ProbMajorMinorOutgroup.txt
        """
