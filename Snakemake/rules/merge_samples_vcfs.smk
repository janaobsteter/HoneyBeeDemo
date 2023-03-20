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
vcfdir = config['vcfDir']

VCFs =[x for x in list_full_paths(vcfdir) if x.endswith(".vcf.gz")]

# creates a list of all files for that chromosome
files_chr = [name in VCFs if '{chromosome}' in name]
noFiles_chr = len(files_chr)
identifier = range(noFiles_chr)

if noFiles_chr == 1:
    print('Nothing to merge')
if noFiles_chr == 0:
    print('Missing chromosome file')
if noFiles_chr > 1:
    print('Merging files...')

    # Get sample list for each vcf
    rule extract_vcf_samples:
        input:
            expand(file, file=files_chr)
        output:
            ids=[expand(vcfdir + 'ids' + i + '.txt' for i in identifier)]
        shell:
            """
            bcftools query -l {input} > {output}
            """

    # Compare sample lists and make sure tehy don't overlap -
    #   keep only firt element if present more than once
    rule remove_duplicates:
        input:
            rule.extract_vcf_samples.output.ids
        output:
            expand(vcfdir + 'ids' + i + 'filtered.txt' for i in identifier)
        run:
            import itertools
            import numpy as np

            for a, b in itertools.combinations(list(input), 2):
                for element in a:
                    if element in b:
                        np.savetxt(output, b.remove(element))

    # Filter vcf based on sample list
    rule filter_samples:
        input:
            vcfs=expand(file, file=files_chr)
            ids=rules.extract_vcf_samples.output
        output:
            vcf=expand(vcfdir + 'Chr{chromosome}_' + i + 'vcf.gz' for i in identifier)
        shell:
            """
            bcftools view -S {input.ids} --force-samples {input.vcfs} -O z {output.vcf}
            bcftools index {output.vcf}
            rm input.ids
            mkdir rawVcfs
            mv {input.vcfs} rawVcfs
            """

    rule bcftools_merge:
        input:
            vcfs=expand(files, files=rules.filter_samples.output.vcf)
            index=expand(indexes, indexes=rules.filter_samples.output.vcf + '.tbi')
        output:
            vcf=[vcfdir + x for x in expand("Chr{{chromosome}}_raw.vcf.gz")]
        shell:
            """
            bcftools bcftools merge -m all {input.vcfs} -O z {output.vcf}
            bcftools index {output.vcf}
            mkdir rawVcfs
            mv {input.vcfs} {input.index} rawVcfs
            """

    rule shapeit_phase:
        input:
            vcf=rules.merge.output.vcf
            map=config['genMaps'] + 'chr{chromosome}.txt'
        output:
            vcf=[vcfdir + x for x in expand("Chr{{chromosome}}.vcf.gz")]
        shell:
            """
            source activate shapeit4am

            shaipeit --input {input.vcf} --map {input.map} --region {wildcards.chromosome} --output {output}
            """
