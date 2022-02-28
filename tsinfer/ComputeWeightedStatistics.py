import pandas as pd
import sys
import json
import tskit
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats
from scipy import cluster
import scipy

filePath = sys.argv[1]
sampleFile = filePath + "/Chr16.trees"
metaFile = sys.argv[2]

meta = pd.read_csv(metaFile)


drone_ts = tskit.load(sampleFile)

sample_nodes = [drone_ts.node(n) for n in drone_ts.samples()]

sample_ids = [n.id for n in sample_nodes]

sample_names = [
    json.loads(drone_ts.individual(n.individual).metadata)["name"]
    for n in sample_nodes
]

sample_pops = [
    json.loads(drone_ts.population(n.population).metadata)["lineage"]
    for n in sample_nodes
]

sample_subspecie = [
    json.loads(drone_ts.population(n.population).metadata)["subspecies"]
    for n in sample_nodes
]

samples_listed_by_population = [
    drone_ts.samples(population=pop_id)
    for pop_id in range(drone_ts.num_populations)
]

samples_listed_by_all = [ind.nodes for ind in drone_ts.individuals()]

sample_ids_names = [
    (n.id, json.loads(drone_ts.individual(n.individual).metadata)["name"])
    for n in sample_nodes
]

sample_id_country = [
    (id, list(meta.Country[meta.ID == name])[0]) if name in list(meta.ID) else (id, "FRAReud") for id, name in sample_ids_names
]

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


def get_lineage_colours():
    return {
        "A": sns.color_palette("Greens", 2)[1],
        "C": sns.color_palette("Blues", 1)[0],
        "H": sns.color_palette("Greys", 3)[0],
        "M": sns.color_palette("Reds", 2)[1],
        "O": sns.color_palette("Purples", 2)[1],
    }

def get_subspecies_colours():
    return {
        "capensis": sns.color_palette("Greens", 2)[1],
        "scutellata": sns.color_palette("Greens", 2)[1],
        "unicolor": sns.color_palette("Greens", 2)[1],
        "ligustica": sns.color_palette("Blues", 1)[0],
        "carnica": sns.color_palette("Blues", 1)[0],
        "hybrid": sns.color_palette("Greys", 3)[0],
        "mellifera": sns.color_palette("Reds", 2)[1],
        "caucasica": sns.color_palette("Purples", 2)[1],
    }



def get_ind_colours(dataframe):
    # TODO add option to give shades for the different pops.
    lineage_colours = get_lineage_colours()
    ind_colour_map = {}
    for lineage, ind in zip(dataframe.population, dataframe.individual):
        ind_colour_map[ind] = lineage_colours[lineage]
    return ind_colour_map

##############################
# Fst
for chromosome in range(1, 17):
    print(chromosome)
    drone_ts = tskit.load(filePath + "/Chr" + str(chromosome) + "_processed.trees")

    lineagePairs = [(0, 1), (0, 2), (0, 3), (0,4),
                    (1, 2), (1, 3), (1, 4),
                    (2,3), (2,4),
                    (3,4)]
    fstLineage = drone_ts.Fst([lineages[lineage] for lineage in lineages.keys()],
                              indexes=lineagePairs)
    fstLineageDF = pd.DataFrame({"Fst": fstLineage})
    fstLineageDF['Subspecie1'] = [list(lineages.keys())[x] for x,y in lineagePairs]
    fstLineageDF['Subspecie2'] = [list(lineages.keys())[y] for x,y in lineagePairs]
    fstLineageDF.to_csv("FstLineages" + str(chromosome) + ".csv", index = False)


    if chromosome == 1:
        fstLineageDFChrPer = fstLineageDF
        fstLineageDFChrPer['Fst'] = fstLineageDFChrPer['Fst'] * chrPer[chromosome - 1]
    else:
        fstLineageDFChrPer['Fst'] = fstLineageDFChrPer['Fst']  + fstLineageDF['Fst'] *  chrPer[chromosome - 1]


    subspeciePairs = [(0, 1), (0, 2), (0, 3), (0,4), (0,5), (0,6), (0,7),
                      (1, 2), (1, 3), (1, 4), (1,5), (1,6), (1,7),
                      (2,3), (2,4), (2,5), (2,6), (2,7),
                      (3,4), (3,5), (3,6), (3,7),
                      (4,5), (4,6), (4,7),
                      (5,6), (5,7),
                      (6,7)]
    fstSubspecie = drone_ts.Fst([subspecies[subspecie] for subspecie in subspecies.keys()],
                                indexes=subspeciePairs)
    fstSubspecieDF = pd.DataFrame({"Fst": fstSubspecie})
    fstSubspecieDF['Subspecie1'] = [list(subspecies.keys())[x] for x,y in subspeciePairs]
    fstSubspecieDF['Subspecie2'] = [list(subspecies.keys())[y] for x,y in subspeciePairs]
    fstSubspecieDF.to_csv("FstSubspecie" + str(chromosome) + ".csv", index = False)


    if chromosome == 1:
        fstSubspecieDFChrPer = fstSubspecieDF
        fstSubspecieDFChrPer['Fst'] = fstSubspecieDFChrPer['Fst'] * chrPer[chromosome - 1]
    else:
        fstSubspecieDFChrPer['Fst'] = fstSubspecieDFChrPer['Fst']  + fstSubspecieDF['Fst'] *  chrPer[chromosome - 1]


    tajimaLineage = drone_ts.Tajimas_D(list(lineages.values()))
    tajimaLineageDF = pd.DataFrame({"TajimaD": tajimaLineage})
    tajimaLineageDF['Lineage'] = lineages.keys()
    tajimaLineageDF.to_csv("TajimaDLineages" + str(chromosome) + ".csv", index = False)


    if chromosome == 1:
        tajimaLineageDFChrPer = tajimaLineageDF
        tajimaLineageDFChrPer['TajimaD'] = tajimaLineageDFChrPer['TajimaD'] * chrPer[chromosome - 1]
    else:
        tajimaLineageDFChrPer['TajimaD'] = tajimaLineageDFChrPer['TajimaD']  + tajimaLineageDF['TajimaD'] *  chrPer[chromosome - 1]


    tajimaLineageDFChrPer.to_csv("TajimaDLineageChrPer.csv", index = False)


    tajimaSubspecie = drone_ts.Tajimas_D(list(subspecies.values()))
    tajimaSubspecieDF = pd.DataFrame({"TajimaD": tajimaSubspecie})
    tajimaSubspecieDF['Subspecie'] = subspecies.keys()


    if chromosome == 1:
        tajimaSubspecieDFChrPer = tajimaSubspecieDF
        tajimaSubspecieDFChrPer['TajimaD'] = tajimaSubspecieDFChrPer['TajimaD'] * chrPer[chromosome - 1]
    else:
        tajimaSubspecieDFChrPer['TajimaD'] = tajimaSubspecieDFChrPer['TajimaD']  + tajimaSubspecieDF['TajimaD'] *  chrPer[chromosome - 1]


    diversityLineage = drone_ts.diversity(list(lineages.values()))
    diversityLineageDF = pd.DataFrame({"Diversity": diversityLineage})
    diversityLineageDF['Lineage'] = lineages.keys()
    diversityLineageDF.to_csv("DiversityLineages" + str(chromosome) + ".csv", index = False)


    if chromosome == 1:
        diversityLineageDFChrPer = diversityLineageDF
        diversityLineageDFChrPer['Diversity'] = diversityLineageDFChrPer['Diversity'] * chrPer[chromosome - 1]
    else:
        diversityLineageDFChrPer['Diversity'] = diversityLineageDFChrPer['Diversity']  + diversityLineageDF['Diversity'] *  chrPer[chromosome - 1]


    diversityLineageDFChrPer.to_csv("DiversityLineageChrPer.csv", index = False)


    diversitySubspecie = drone_ts.diversity(list(subspecies.values()))
    diversitySubspecieDF = pd.DataFrame({"Diversity": diversitySubspecie})
    diversitySubspecieDF['Subspecie'] = subspecies.keys()


    if chromosome == 1:
        diversitySubspecieDFFChrPer = diversitySubspecieDF
        diversitySubspecieDFFChrPer['Diversity'] = diversitySubspecieDFFChrPer['Diversity'] * chrPer[chromosome - 1]
    else:
        diversitySubspecieDFFChrPer['Diversity'] = diversitySubspecieDFFChrPer['Diversity']  + diversitySubspecieDF['Diversity'] *  chrPer[chromosome - 1]


    gnn = drone_ts.genealogical_nearest_neighbours(
        drone_ts.samples(), samples_listed_by_all
    )


    cols = {drone_ts.node(u).individual: gnn[:, u] for u in drone_ts.samples()}
    cols["population"] = [json.loads(drone_ts.population(drone_ts.node(u).population).metadata)["lineage"] for u in drone_ts.samples()]
    cols["subspecie"] = [json.loads(drone_ts.population(drone_ts.node(u).population).metadata)["subspecies"] for u in drone_ts.samples()]
    cols["individual"] = [drone_ts.node(u).individual for u in drone_ts.samples()]
    df = pd.DataFrame(cols)
    df.to_csv(filePath + "Statistics/GnnInd_Chr" + str(chromosome) + ".csv")

    chrGnn = pd.read_csv(filePath +'Statistics/GnnInd_Chr' + str(chromosome) + '.csv', index_col = None)

    if chromosome == 1:
        combinedGnn = chrGnn.iloc[:, 1:692]
    else:
        combinedGnn = combinedGnn + chrGnn.iloc[:, 1:692] * chrPer[chromosome - 1]

fstLineageDFChrPer.to_csv(filePath +"Statistics/FstLineageChrPer.csv", index = False)
fstSubspecieDFChrPer.to_csv(filePath +"Statistics/FstSubspecieChrPer.csv", index = False)
tajimaSubspecieDF.to_csv(filePath +"Statistics/TajimaDSubspecie" + str(chromosome) + ".csv", index = False)
tajimaSubspecieDFChrPer.to_csv(filePath +"Statistics/TajimaDSubspecieChrPer.csv", index = False)
diversitySubspecieDF.to_csv(filePath +"Statistics/DiversityDSubspecie" + str(chromosome) + ".csv", index = False)
diversitySubspecieDFFChrPer.to_csv(filePath +"Statistics/DiversitySubspecieChrPer.csv", index = False)

combinedGnn.to_csv(filePath +"Statistics/WeightedGnn_Ind.csv", index = False)

#combinedGnn = pd.read_csv(filePath +"Statistics/WeightedGnn_Ind.csv")

# Plot the data
figsize = (10, 10)
combinedGnnG = combinedGnn.copy() #.groupby("individual").mean()
combinedGnn['subspecie'] = sample_subspecie
combinedGnn['individual'] = chrGnn.individual
combinedGnn['population'] = chrGnn.population
combinedGnn["country"] = [country for (id, country) in sample_id_country]
# logit normalise
for col in list(combinedGnnG):
    combinedGnnG[col] = [scipy.special.logit(x+ 0.0001) for x in list(combinedGnnG[col])]
row_linkage = scipy.cluster.hierarchy.linkage(combinedGnnG, method="average")
order = scipy.cluster.hierarchy.leaves_list(row_linkage)
x_pop = combinedGnnG.index.values[order]
combinedGnnG.columns = [int(x) for x in combinedGnnG.columns]

#Transform
cg = sns.clustermap(
    combinedGnnG[x_pop], row_linkage=row_linkage, col_cluster=False, # dfg[x_pop] dfgOrder[list(dfgOrder.individual)]
    figsize=figsize, rasterized=True, row_colors = pd.Series(get_ind_colours(combinedGnn)))
cg.ax_heatmap.set_yticks([])
cg.ax_heatmap.set_xticks([])
#cg.cax.set_visible(False)

from matplotlib.patches import Patch

lineages = ['A', 'M', 'C', 'O', 'H']
lineages = get_lineage_colours().keys()
lineagesCols = get_lineage_colours()
handles = [Patch(facecolor=lineagesCols[name]) for name in lineages]
plt.legend(handles, lineagesCols, title='Subspecies',
           bbox_to_anchor=(1, 1), bbox_transform=plt.gcf().transFigure, loc='upper right',
           prop={'size': 20}, ncol=3, title_fontsize=20)
subspecies = get_subspecies_colours().keys()
subspeciesCols = get_subspecies_colours()
handles = [Patch(facecolor=subspeciesCols[name]) for name in subspecies]
plt.legend(handles, subspeciesCols, title='Subspecies',
           bbox_to_anchor=(1, 1), bbox_transform=plt.gcf().transFigure, loc='upper right',
           prop={'size': 20}, ncol=3, title_fontsize=20)

plt.show()
cg.figure.savefig(filePath + "Statistics/GnnInd_Logit_Lineages.png")
cg.figure.savefig(filePath + "Statistics/GnnInd_Logit_SubSpecies.png")
