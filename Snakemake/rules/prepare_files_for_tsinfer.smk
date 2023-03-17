from os import listdir
configfile: "../config/tsinfer.yaml"

rule all:
    input:
        expand("Tsinfer/Chr{chromosome}.vcf.gz", chr = range(1, config['noChromosomes'] + 1))

noVCFs = length(listdir(config['vcfDir']))
if noVCFs == config['noChromosomes']:
    print("Your VCFs are split.")

if (noVCFs != config['noChromosomes']) & noVCFs != 1:
    print("Check you VCFs.")

if noVCFs == 1:
    rule split_vcfs: #Input VCFs need ot be compressed and indexes
        input:
            vcf:listdir(config['vcfDir'])[0]
            vcfDir: config['vcfDir']
        output:
            vcf:expand("{input.vcfDir}/Chr{chromosome}.vcf.gz", chromosome = range(1, config['noChromosomes'] + 1))
            index:expand("{input.vcfDir}/Chr{chromosome}.vcf.gz.csi", chromosome = range(1, config['noChromosomes'] + 1))
        #envmodules:
        #    config['bcftoolsModule']
        shell:
            """
            bcftools view -r {wildcards.chromosome} -O z {input} > {output}
            bcftools index {output}
            """

rule decompress:
    input:
        "Tsinfer/Chr{chromosome}.vcf.gz" #This will take
    output:
        "Tsinfer/Chr{chromosome}.vcf"
    shell:
        "gunzip {input} {output}"


rule extract_vcf_pos:
    input:
        rules.decompress.output
    output:
        "Tsinfer/VcfPos{chromosome}.txt"
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
        vcfPos=rules.extract_vcf_pos.output# This has more lines
        ancestralAllele=config['ancestralAllele'] # The ancestral file has to have chr_pos and AA, split with a tab
    output:
        "Tsinfer/AncestralVcfMatch{chromosome}.txt" # THis files needs to contain all the variants from the vcf (blank space)
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
        vcf=rules.decompress.output
        ancestralAllele=rules.match_ancestral_vcf.output
    output:
        "Tsinfer/Chr{chromosome}_ancestral.vcf"
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
        rules.change_infoAA_vcf.output
    output:
        "Tsinfer/Chr{chromosome}_ancestral.vcf.gz"
    shell:
        "bgzip {input}"

rule index_vcf:
    input:
        rules.compress_vcf.output
    output:
        "Tsinfer/Chr{chromosome}_ancestral.vcf.gz.csi"
    shell:
        "bcftools index -f {input} > {output}"
