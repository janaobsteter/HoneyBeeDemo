#
# Files must be at the vcf folder compressed and indexed
#
import os
configfile: "../config/tsinfer.yaml"

def list_full_paths(directory):
    return [os.path.join(directory, file) for file in os.listdir(directory)]

rule all:
    input:
        expand(contig['vcfDir'] + "Chr{chromosome}.vcf.gz", chromosome = range(1, config['noChromosomes'] + 1))

workdir: config['workdir']

VCFs =[x for x in list_full_paths(config['vcfDir']) if x.endswith(".vcf.gz")]

# creates a list of all files for that chromosome
files_chr = [name in VCFs if '{chromosome}' in name]

if len(files_chr) == 1:
    print('Nothing to merge')
if len(files_chr) == 0:
    print('Missing chromosome file')
if len(files_chr) > 1:
    print('Merging files...')

    rule extract_vcf_samples:
        input:
            expand(files_chr)
        output:
            contig['vcfDir'] + ids.txt
        shell:
            """
            bcftools query -l {input} >> tmp
            sort tmp | uniq -u > {output}
            """

    rule filter_samples:
        input:
            vcfs=expand(files_chr)
            ids=rules.extract_vcf_samples.output
        output:
            vcf=expand(contig['vcfDir'] + 'Chr{chromosome}_{n}.vcf.gz', n=range(len(files_chr)))
        shell:
            """
            bcftools view -S {input.ids} --force-samples {input.vcfs} -O z {output.vcf}
            bcftools index {output.vcf}
            rm input.ids
            mkdir rawVcfs
            mv {input.vcfs} rawVcfs
            """

    rule merge:
        input:
            vcfs=expand(rules.filter_samples.output.vcf)
            index=expand(rules.filter_samples.output.vcf + '.tbi')
        output:
            vcf=contig['vcfDir'] + 'Chr{chromosome}.vcf.gz'
        shell:
            """
            bcftools bcftools merge -m all {input.vcfs} -O z {output.vcf}
            bcftools index {output.vcf}
            mkdir rawVcfs
            mv {input.vcfs} {input.index} rawVcfs
            """
