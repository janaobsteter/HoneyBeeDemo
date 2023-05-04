include: ancestral_inference.smk
include: multiple_genome_alignment.smk
include: qc_vcf.smk
include: tsinfer.smk

configfile: "config/tsinfer.yaml"

rule compute_tree_stats:
    input:
        expand("Tsinfer/Chr{chromosome}.trees", chromosome = range(1, config['noChromosomes'] + 1))
    input:
        "Tsinfer/Chr{chromosome}.trees"
    output:
        fst="Tsinfer/Stats/Fst{chromosome}.csv",
        tajimasd="Tsinfer/Stats/TajimasD{chromosome}.csv",
        divergence="Tsinfer/Stats/Divergence{chromosome}.csv",
        diversity="Tsinfer/Stats/Diversity{chromosome}.csv",
        f2="Tsinfer/Stats/F2{chromosome}.csv",
        gRel="Tsinfer/Stats/GeneticRelatednes{chromosome}.csv",
        meanDesc="Tsinfer/Stats/MeanDesc{chromosome}.csv",
        mutations="Tsinfer/Stats/Mutations{chromosome}.csv",
        nodes="Tsinfer/Stats/Nodes{chromosome}.csv",
        sites="Tsinfer/Stats/Sites{chromosome}.csv",
        edges="Tsinfer/Stats/Edges{chromosome}.csv"
    script:
        "scripts/TsinferStats.py"

rule combine_stats:
    input:
        fst=expand("Tsinfer/Stats/Fst{chromosome}.csv", chromosome = range(1:config['noChromosomes'])),
        tajimasd=expand("Tsinfer/Stats/TajimasD{chromosome}.csv", chromosome = range(1:config['noChromosomes'])),
        divergence=expand("Tsinfer/Stats/Divergence{chromosome}.csv", chromosome = range(1:config['noChromosomes'])),
        diversity=expand("Tsinfer/Stats/Diversity{chromosome}.csv", chromosome = range(1:config['noChromosomes'])),
        f2=expand("Tsinfer/Stats/F2{chromosome}.csv", chromosome = range(1:config['noChromosomes'])),
        gRel=expand("Tsinfer/Stats/GeneticRelatednes{chromosome}.csv", chromosome = range(1:config['noChromosomes'])),
        meanDesc=expand("Tsinfer/Stats/MeanDesc{chromosome}.csv", chromosome = range(1:config['noChromosomes'])),
        mutations=expand("Tsinfer/Stats/Mutations{chromosome}.csv", chromosome = range(1:config['noChromosomes'])),
        nodes=expand("Tsinfer/Stats/Nodes{chromosome}.csv", chromosome = range(1:config['noChromosomes'])),
        sites=expand("Tsinfer/Stats/Sites{chromosome}.csv", chromosome = range(1:config['noChromosomes'])),
        edges=expand("Tsinfer/Stats/Edges{chromosome}.csv", chromosome = range(1:config['noChromosomes']))
    output:
        fst="Tsinfer/Stats/FstCombined.csv",
        tajimasd="Tsinfer/Stats/TajimasDCombined.csv",
        divergence="Tsinfer/Stats/DivergenceCombined.csv",
        diversity="Tsinfer/Stats/DiversityCombined.csv",
        f2="Tsinfer/Stats/F2Combined.csv",
        gRel="Tsinfer/Stats/GeneticRelatednesCombined.csv",
        meanDesc="Tsinfer/Stats/MeanDescCombined.csv",
        mutations="Tsinfer/Stats/MutationsCombined.csv",
        nodes="Tsinfer/Stats/NodesCombined.csv",
        sites="Tsinfer/Stats/SitesCombined.csv",
        edges="Tsinfer/Stats/EdgesCombined.csv"
    script:
        "scripts/CombineStats.py"
