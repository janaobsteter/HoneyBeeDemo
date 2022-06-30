include: multiple_genome_alignment.smk
include: qc_vcf.smk

configfile: "config/ancestral.yaml"

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
        alignedAlleles:"AlignedSnps_focal_CarnicaChrSorted.txt"
        vcf:config['rawVcf']
        out: "VcfInfo"
    output:
        "VcfInfo.INFO"
    script:
        "scripts/ExtractPosVcf.sh {input.vcf} {input.alignedPos} {input.out}"

rule extract_pos_from_aligned_alleles_from_vcf:
    input:
        "VcfInfo.INFO"
    output:
        "VcfFullPos.txt"
    shell:
        "awk '{print $1 "_" $2}' {input} > {output}"

rule extract_vcf_alleles_from_aligned:
    input:
        alignedAlleles:"AlignedSnps_focal_CarnicaChrSorted.txt"
        vcfPos:"VcfFullPos.txt"
    output:
        "AlignedSnps_focal_CarnicaChrSortedVcf.txt"
    shell:
        "grep -Fwf {input.vcfPos} {input.alignedPop} > {output}"

rule create_estsfs_dicts:
    input:
        alignedAlleles: "AlignedSnps_focal_CarnicaChrSortedVcf.txt"
        vcfAlleles: "VcfInfo.INFO"
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
    shell:
        """
        for cycle in $(seq {config['noEstSfsChunks']})
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

rule run_estsfs:
    input:
        config="config-kimura_3o.txt",
        dict="EstSfs_Dict{chunk}E.csv",
        seed="seedfile.txt"
    output:
        text="EstsfsOutput/outputEtsfs{chunk}.txt",
        pvalue="EstsfsOutput/output-{chunk}-pvalues.txt"
    shell:
        """
        mkdir EstsfsOutput
        ./estsfs {input.config} {input.dict} {output.text} {output.pvalue}
        """

rule determine_ancestral_allele:
    input:
        "EstsfsOutput/output-{chunk}-pvalues.txt"
    output:
        "EstsfsOutput/AncestralAllele{chunk}_3o.csv"
    script:
        "scripts/DetermineAncestralAllele.py"

rule combine_estsfs_dicts:
    input:
        expand("EstSfs_Dict{chunk}.csv", chunk = [x for x in range(config['noEstSfsChunks'])])
