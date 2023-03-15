include: ancestral_inference.smk
include: multiple_genome_alignment.smk
include: qc_vcf.smk

configfile: "config/tsinfer.yaml"

rule all:
    input:
        expand("Tsinfer/Chr{chromosome}.trees", chromosome = range(1, config['noChromosomes'] + 1)

rule extract_vcf_pos:
    input:
        config['vcf']
    output:
        "Tsinfer/vcfPos.txt"
    #conda:
    #    config['tsinferEnv']
    envmodules:
        config['bcftoolsModule']
    shell:
    "bcftools query -f '%CHROM %POS\n' {input} > {output}"

rule match_ancestral_vcf:
    input:
        vcfPos="Tsinfer/vcfPos.txt"# This has more lines
        ancestralAllele=config['ancestralAllele']
    output:
        "Tsinfer/AncestralVcfMatch.txt"
    shell:
    """
    for line in $(cat {input.vcf}});
    do
      grep $line {input.ancestralAllele} || echo "";
    done
    """

rule change_infoAA_vcf:
    input:
        vcf=config['vcf']
        ancestralAllele="Tsinfer/AncestralVcfMatch.txt"
    output:
        "Tsinfer/Vcf_AncestralInfo.vcf"
    shell:
        """
        NUM=$(( $(grep "##" {input.vcf} | wc -l) + 1 ))
        awk -v NUM=$NUM 'NR==FNR{a[NR] = $2; next} FNR<29{print}; \
        FNR==5{printf "##INFO=<ID=AA,Number=1,Type=String,Description=\"Ancestral Allele\">\n"}; \
        FNR>28{$8=a[FNR-28]; print}' OFS="\t" All_AAInfo.txt test.vcf awk '{FNR>28{print a[FNR-28]}}' \
        {input.ancestralAllele} {input.vcf} > {output}
        """

rule split_vcf:
    input:
        "Tsinfer/Vcf_AncestralInfo.vcf"
    output:
        expand("Tsinfer/Chr{{chromosome}}.vcf")
    shell:
        "bcftools view -r {wildcards.chromosome} {input} > {output}"

rule prepare_sample_file:
    input:
        vcf="Tsinfer/Chr{chromosome}.vcf"
        meta=config['meta']
    output:
        "Tsinfer/Chr{chromosome}.samples"
    script:
        "scripts/PrepareTsinferSampleFile.py"

rule infer:
    input:
        "Tsinfer/Chr{chromosome}.samples"
    output:
        "Tsinfer/Chr{chromosome}.trees"
    script:
        "scripts/InferTrees.py"
