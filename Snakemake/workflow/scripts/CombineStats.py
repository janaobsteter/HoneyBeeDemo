import pandas as pd
import os

fstFiles = [x for x in snakemake.input['fst']]
tajimasdFiles = [x for x in snakemake.input['tajimasd']]
divergenceFiles = [x for x in snakemake.input['divergence']]
diversityFiles = [x for x in snakemake.input['diversity']]
f2Files = [x for x in snakemake.input['f2']]
gRelFiles = [x for x in snakemake.input['gRel']]
meanDescFiles = [x for x in snakemake.input['meanDesc']]
mutationsFiles = [x for x in snakemake.input['mutations']]
nodesFiles = [x for x in snakemake.input['nodes']]
sitesFiles = [x for x in snakemake.input['sites']]
edgesFiles = [x for x in snakemake.input['edges']]

genome = {
    "assembly_accession": "GCA_003254395.2",
    "assembly_name": "Amel_HAv3.1",
    "chromosomes": {
        "CM009931.2": {"length": 27754200, "synonyms": ["NC_037638.1"]},
        "CM009932.2": {"length": 16089512, "synonyms": ["NC_037639.1"]},
        "CM009933.2": {"length": 13619445, "synonyms": ["NC_037640.1"]},
        "CM009934.2": {"length": 13404451, "synonyms": ["NC_037641.1"]},
        "CM009935.2": {"length": 13896941, "synonyms": ["NC_037642.1"]},
        "CM009936.2": {"length": 17789102, "synonyms": ["NC_037643.1"]},
        "CM009937.2": {"length": 14198698, "synonyms": ["NC_037644.1"]},
        "CM009938.2": {"length": 12717210, "synonyms": ["NC_037645.1"]},
        "CM009939.2": {"length": 12354651, "synonyms": ["NC_037646.1"]},
        "CM009940.2": {"length": 12360052, "synonyms": ["NC_037647.1"]},
        "CM009941.2": {"length": 16352600, "synonyms": ["NC_037648.1"]},
        "CM009942.2": {"length": 11514234, "synonyms": ["NC_037649.1"]},
        "CM009943.2": {"length": 11279722, "synonyms": ["NC_037650.1"]},
        "CM009944.2": {"length": 10670842, "synonyms": ["NC_037651.1"]},
        "CM009945.2": {"length": 9534514, "synonyms": ["NC_037652.1"]},
        "CM009946.2": {"length": 7238532, "synonyms": ["NC_037653.1"]},
        "CM009947.2": {"length": 16343, "synonyms": ["NC_001566.1", "MT"]},
    },
}

chromosomes = genome['chromosomes']
chrLengths = [value['length'] for value in chromosomes.values()]
chrPer = [x / sum(chrLengths) for x in chrLengths]

for chromosome in range(snakemake.config['noChromosomes']):
    fstDF = pd.read_csv(stFiles[chromosome])
    tajimasDF = pd.read_csv(tajimasdFiles[chromosome])
    divergenceDF = pd.read_csv(divergenceFiles[chromosome])
    diversityDF = pd.read_csv(diversityFiles[chromosome])
    f2DF = pd.read_csv(f2Files[chromosome])
    gRelDF = pd.read_csv(gRelFiles[chromosome])
    meanDescDF = pd.read_csv(meanDescFiles[chromosome])
    if chromosome == 0:
        #Fst
        fstCombined = fstDF
        fstCombined['Fst'] = fstCombined['Fst'] * chrPer[chromosome]
        #TajimasD
        tajimasDCombined = tajimasDF
        tajimasDCombined['TajimasD'] = tajimasDCombined['TajimasD'] * chrPer[chromosome]
        #Divergence
        divergenceCombined = divergenceDF
        divergenceCombined['Divergence'] = divergenceCombined['Divergence'] * chrPer[chromosome]
        #Diversity
        diversityCombined = divergenceDF
        divergenceCombined['Diversity'] = divergenceCombined['Diversity'] * chrPer[chromosome]
        #Diversity
        diversityCombined = diversityDF
        diversityCombined['Diversity'] = diversityCombined['Diversity'] * chrPer[chromosome]
        #f2
        f2Combined = f2DF
        f2Combined['f2'] = f2Combined['f2'] * chrPer[chromosome]
        #Genetic relatedness
        gRelCombined = gRelDF
        gRelCombined['GeneticRel'] = gRelCombined['GeneticRel'] * chrPer[chromosome]
        #Mean descendants
        meanDescDF.loc[:, "Chromosome" = str(chromosome + 1)]
        meanDescCombined = meanDescDF
    else:
        fstCombined['Fst'] = fstCombined['Fst']  + fstDF['Fst'] *  chrPer[chromosome]
        tajimasDCombined['TajimasD'] = tajimasDCombined['TajimasD']  + tajimasDF['TajimasD'] *  chrPer[chromosome]
        divergenceCombined['Divergence'] = divergenceCombined['Divergence']  + divergenceDF['Divergence'] *  chrPer[chromosome]
        diversityCombined['Diversity'] = diversityCombined['Diversity']  + diversityDF['Diversity'] *  chrPer[chromosome]
        f2Combined['f2'] = f2Combined['f2']  + f2DF['f2'] *  chrPer[chromosome]
        gRelCombined['GeneticRel'] = gRelCombined['GeneticRel']  + gRelDF['GeneticRel'] *  chrPer[chromosome]
        meanDescDF.loc[:, "Chromosome" = str(chromosome + 1)]
        meanDescCombined = pd.concat([meanDescCombined, meanDescDF])


fstCombined.to_csv(snakemake.output['fst'], index = False)
tajimasDCombined.to_csv(snakemake.output['tajimasd'], index = False)
divergenceCombined.to_csv(snakemake.output['divergence'], index = False)
diversityCombined.to_csv(snakemake.output['diversity'], index = False)
f2Combined.to_csv(snakemake.output['f2'], index = False)
gRelCombined.to_csv(snakemake.output['gRel'], index = False)
meanDescCombined.to_csv(snakemake.output['meanDesc'], index = False)

# Mutations, edges,sites, nodes
os.system("{ head -n1 " + mutationFiles[0] + "; for f in " + " ".join(mutationFiles) + """; do tail -n+2 "$f"; done; } > """ + snakemake.output['mutations'])
os.system("{ head -n1 " + edgesFiles[0] + "; for f in " + " ".join(edgesFiles) + """; do tail -n+2 "$f"; done; } > """ + snakemake.output['edges'])
os.system("{ head -n1 " + nodesFiles[0] + "; for f in " + " ".join(nodesFiles) + """; do tail -n+2 "$f"; done; } > """ + snakemake.output['nodes'])
os.system("{ head -n1 " + sitesFiles[0] + "; for f in " + " ".join(sitesFiles) + """; do tail -n+2 "$f"; done; } > """ + snakemake.output['sites'])
