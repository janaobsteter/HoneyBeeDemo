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
        noChunks: "1000"
        alignedAlleles: "AlignedSnps_focal_CarnicaChrSortedVcf.txt"
        vcfAlleles: "VcfInfo.INFO"
    output:
        expand("EstSfs_Dict{chunk}.csv", chunk = [x for x in range(noChunks)])
    script:
        "scripts/CreateInputForEstsfs_fromWGAbed_Cactus.py"
        #CreateInputForEstsfs_Loop.sh This is qsub

rule edit_estsfs_dicts:
    input:
        expand("EstSfs_Dict{chunk}.csv", chunk = [x for x in range(noChunks)])
