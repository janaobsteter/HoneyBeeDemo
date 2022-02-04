import tskit
import json
import pandas as pd
import seaborn as sns
from collections import Counter



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


drone_ts = tskit.load("Trees/Chr16_processed.trees")

# Obtain node info - but just for samples!
# There is also intermediate nodes!!! (total #nodes > #nodes for samples)
sample_nodes = [drone_ts.node(n) for n in drone_ts.samples()]

# Get samples ids
sample_ids = [n.id for n in sample_nodes]
sample_nodes = [drone_ts.node(n) for n in drone_ts.samples()]
# Get sample names
sample_names = [
    json.loads(drone_ts.individual(n.individual).metadata)["name"]
    for n in sample_nodes
]

# Get sample population
sample_pops = [
    json.loads(drone_ts.population(n.population).metadata)["lineage"]
    for n in sample_nodes
]
Counter(sample_pops)

# Get sample population
sample_subspecie = [
    json.loads(drone_ts.population(n.population).metadata)["subspecies"]
    for n in sample_nodes
]
Counter(sample_subspecie)

samples_listed_by_population = [
    drone_ts.samples(population=pop_id)
    for pop_id in range(drone_ts.num_populations)
]

samples_listed_by_all = [ind.nodes for ind in drone_ts.individuals()]

lineages = dict()
for (node, pop) in zip(sample_ids, sample_pops):
    if not pop in lineages.keys():
        lineages[pop] = [node]
    else:
        lineages[pop].append(node)

subspecies = dict()
for (node, pop) in zip(sample_ids, sample_subspecie):
    if not pop in subspecies.keys():
        subspecies[pop] = [node]
    else:
        subspecies[pop].append(node)


for chromosome in range(1, 17):
    print(chromosome)
    drone_ts = tskit.load("Trees/Chr" + str(chromosome) + "_processed.trees")

    lineagePairs = [(0, 1), (0, 2), (0, 3), (0,4),
                      (1, 2), (1, 3), (1, 4),
                      (2,3), (2,4),
                      (3,4)]
    f2Lineage = drone_ts.f2([lineages[lineage] for lineage in lineages.keys()],
                   indexes=lineagePairs)
    f2LineageDF = pd.DataFrame({"f2": f2Lineage})
    f2LineageDF['Subspecie1'] = [list(lineages.keys())[x] for x,y in lineagePairs]
    f2LineageDF['Subspecie2'] = [list(lineages.keys())[y] for x,y in lineagePairs]
    f2LineageDF.to_csv("F2Lineages" + str(chromosome) + ".csv", index = False)


    if chromosome == 1:
        f2LineageDFChrPer = f2LineageDF
        f2LineageDFChrPer['f2'] = f2LineageDFChrPer['f2'] * chrPer[chromosome - 1]
    else:
        f2LineageDFChrPer['f2'] = f2LineageDFChrPer['f2']  + f2LineageDF['f2'] *  chrPer[chromosome - 1]


    f2LineageDFChrPer.to_csv("F2LineageChrPer.csv", index = False)


    subspeciePairs = [(0, 1), (0, 2), (0, 3), (0,4), (0,5), (0,6), (0,7),
                      (1, 2), (1, 3), (1, 4), (1,5), (1,6), (1,7),
                      (2,3), (2,4), (2,5), (2,6), (2,7),
                      (3,4), (3,5), (3,6), (3,7),
                      (4,5), (4,6), (4,7),
                      (5,6), (5,7),
                      (6,7)]
    f2Subspecie = drone_ts.f2([subspecies[subspecie] for subspecie in subspecies.keys()],
                              indexes=subspeciePairs)
    f2SubspecieDF = pd.DataFrame({"f2": f2Subspecie})
    f2SubspecieDF['Subspecie1'] = [list(subspecies.keys())[x] for x,y in subspeciePairs]
    f2SubspecieDF['Subspecie2'] = [list(subspecies.keys())[y] for x,y in subspeciePairs]
    f2SubspecieDF.to_csv("F2Subspecie" + str(chromosome) + ".csv", index = False)


    if chromosome == 1:
        f2SubspecieDFChrPer = f2SubspecieDF
        f2SubspecieDFChrPer['f2'] = f2SubspecieDFChrPer['f2'] * chrPer[chromosome - 1]
    else:
        f2SubspecieDFChrPer['f2'] = f2SubspecieDFChrPer['f2']  + f2SubspecieDF['f2'] *  chrPer[chromosome - 1]

    f2SubspecieDFChrPer.to_csv("F2SubspecieChrPer.csv", index = False)

