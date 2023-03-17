#include: ancestral_inference.smk
#include: multiple_genome_alignment.smk
#include: qc_vcf.smk

configfile: "../config/tsinfer.yaml"
workdir: config['workdir']

rule all:
    input:
        expand("Tsinfer/Chr{chromosome}.trees", chromosome = range(1, config['noChromosomes'] + 1))

rule extract_vcf_pos:
    input:
        config['vcf']
    output:
        "Tsinfer/VcfPos.txt"
    #conda:
    #    config['tsinferEnv']
    #envmodules:
    #    config['bcftoolsModule']
    shell:
        """
        bcftools query -f '%CHROM %POS\n' {input} > tmp
        awk '{{print $1"_"$2}}' tmp > {output}
        """

rule match_ancestral_vcf:
    input:
        vcfPos="Tsinfer/VcfPos.txt",# This has more lines
        ancestralAllele=config['ancestralAllele'] # The ancestral file has to have chr_pos and AA, split with a tab
    output:
        "Tsinfer/AncestralVcfMatch.txt"
    shell:
        """
        for line in $(cat {input.vcfPos});
        do
          grep $line {input.ancestralAllele} || echo "";
        done > tmp
        awk -F"," '{{print $1"\t"$2}}' tmp > {output}
        """

rule change_infoAA_vcf:
    input:
        vcf=config['vcf'],
        ancestralAllele="Tsinfer/AncestralVcfMatch.txt"
    output:
        "Tsinfer/Vcf_AncestralInfo.vcf"
    shell:
        """
        HEADERNUM=$(( $(grep "##" {input.vcf} | wc -l) + 1 ))
        INFOLINE=$(( $(grep -Fn "INFO" {input.vcf} | cut --delimiter=":" --fields=1  | head -n1) ))
        awk -v HEADER=$HEADERNUM -v INFO=$INFOLINE 'NR==FNR{{{{a[FNR] = $2; next}}}} FNR<=HEADER{{{{print}}}}; \
        FNR==INFO{{{{printf "##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele">\\n"}}}}; \
        FNR>HEADER{{{{$8=a[FNR-HEADER]; print}}}}' OFS="\t" {input.ancestralAllele} {input.vcf} > {output}
        """
rule compress_vcf:
    input:
        "Tsinfer/Vcf_AncestralInfo.vcf"
    output:
        "Tsinfer/Vcf_AncestralInfo.vcf.gz"
    shell:
        "bgzip {input}"

rule index_vcf:
    input:
        "Tsinfer/Vcf_AncestralInfo.vcf.gz"
    output:
        "Tsinfer/Vcf_AncestralInfo.vcf.gz.csi"
    shell:
        "bcftools index -f {input} > {output}"

rule split_vcf:
    input:
        vcf="Tsinfer/Vcf_AncestralInfo.vcf.gz",
        index="Tsinfer/Vcf_AncestralInfo.vcf.gz.csi"
    output:
        expand("Tsinfer/Chr{{chromosome}}.vcf")
    shell:
        "bcftools view -r {wildcards.chromosome} {input} > {output}"

rule prepare_sample_file:
    input:
        vcf="Tsinfer/Chr{chromosome}.vcf",
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
