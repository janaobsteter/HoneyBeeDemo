import tskit
import json
import numpy as np
import pandas as pd
import scipy
import seaborn as sns
from scipy import stats
from scipy import cluster
import seaborn as sns


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



def get_subspecies_colours(dataframe):
    # TODO add option to give shades for the different pops.
    lineage_colours = get_lineage_colours()
    subspecies_colour_map = {}
    for lineage, subspecie in zip(dataframe.Lineage, dataframe.Subspecie):
        subspecies_colour_map[subspecie] = lineage_colours[lineage]
    return subspecies_colour_map


def get_ind_colours(dataframe):
    # TODO add option to give shades for the different pops.
    lineage_colours = get_lineage_colours()
    ind_colour_map = {}
    for lineage, ind in zip(dataframe.population, dataframe.individual):
        ind_colour_map[ind] = lineage_colours[lineage]
    return ind_colour_map



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

# Obtain node info - but just for samples!
# There is also intermediate nodes!!! (total #nodes > #nodes for samples)
drone_ts = tskit.load("Trees/Chr16_processed.trees")

sample_nodes = [drone_ts.node(n) for n in drone_ts.samples()]

# Get samples ids
sample_ids = [n.id for n in sample_nodes]

# Get sample names
sample_names = [
    json.loads(drone_ts.individual(n.individual).metadata)["name"]
    for n in sample_nodes
]


# Get sample population
sample_subspecie = [
    json.loads(drone_ts.population(n.population).metadata)["subspecies"]
    for n in sample_nodes
]

# Get sample population
sample_pops = [
    json.loads(drone_ts.population(n.population).metadata)["lineage"]
    for n in sample_nodes
]



samples_listed_by_population = [
    drone_ts.samples(population=pop_id)
    for pop_id in range(drone_ts.num_populations)
]


samples_listed_by_population_nonH = [
   samples_listed_by_population[x] for x in [0,1,2,3,4,6,7]
]

samples_listed_by_population_H = [
   samples_listed_by_population[x] for x in [5]
][0]


sample_ids_names_H = [
    (n.id, json.loads(drone_ts.individual(n.individual).metadata)["name"])
    for n in sample_nodes if n.id in lineages["H"]
]


samples_listed_by_all_H = [ind.nodes for ind in drone_ts.individuals() if ind.id in lineages["H"]]

meta = pd.read_csv("Meta.csv")
country = [list(meta.Country[meta.ID == name])[0] if name in list(meta.ID) else "NA" for id, name in sample_ids_names_H ]


for chromosome in [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]:
    print(chromosome)
    drone_ts = tskit.load("Trees/Chr" + str(chromosome) + "_processed.trees")

    gnn = pd.DataFrame(drone_ts.genealogical_nearest_neighbours
                       (drone_ts.samples(), samples_listed_by_all_H))

    gnn = gnn.iloc[lineages["H"]]
    gnn.reset_index()
    namesH = [name for (name, ID) in zip(sample_names, sample_ids) if ID in lineages['H']]
    gnn.loc[:, "Names"] = namesH
    lineagesH = [pop for (pop, ID) in zip(sample_pops, sample_ids) if ID in lineages['H']]
    gnn.loc[:, "Lineage"] = lineagesH
    subspeciesH = [pop for (pop, ID) in zip(sample_subspecie, sample_ids) if ID in lineages['H']]
    gnn.loc[:, "Subspecie"] = subspeciesH
    gnn = gnn.drop_duplicates()

    gnn.to_csv("GnnInd_Chr" + str(chromosome) + "_hybrid.csv")



for chromosome in [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]:
    chrGnn = pd.read_csv('GnnInd_Chr' + str(chromosome) + '_hybrid.csv', index_col = None)
    if chromosome == 1:
        combinedGnn = chrGnn.iloc[:, 1:452]
    else:
        combinedGnn = combinedGnn + chrGnn.iloc[:, 1:9] * chrPer[chromosome - 1]

colourPallete = sns.color_palette("Paired", len(country))
colourPallete[14] = (0,0,0)
country_colours = dict(zip(country, colourPallete))
country_plot_colours = [country_colours[country] for (id,name), country in zip(sample_ids_names_H, country)]


# Plot the data
figsize = (10, 10)
combinedGnnG = combinedGnn
combinedGnnG.reset_index()
# Zscore normalise
for col in list(combinedGnnG):
    combinedGnnG[col] = scipy.stats.zscore(combinedGnnG[col])
row_linkage = scipy.cluster.hierarchy.linkage(combinedGnnG, method="average")
order = scipy.cluster.hierarchy.leaves_list(row_linkage)
x_pop = combinedGnnG.index.values[order]

cg = sns.clustermap(
    combinedGnnG, row_linkage=row_linkage, col_cluster=False, # dfg[x_pop] dfgOrder[list(dfgOrder.individual)]
    figsize=figsize, rasterized=True, row_colors = pd.Series(country_plot_colours))
plt.show()





for chromosome in range(1, 17):
    print(chromosome)
    drone_ts = tskit.load("Trees/Chr" + str(chromosome) + "_processed.trees")

    gnn = pd.DataFrame(drone_ts.genealogical_nearest_neighbours
                       (lineages["H"], sample_sets=[lineages["A"], lineages["C"], lineages["M"], lineages["O"]]),
                       columns=["A", "C", "M", "O"])

    namesH = [name for (name, ID) in zip(sample_names, sample_ids) if ID in lineages['H']]

    gnn.loc[:, "Names"] = namesH
    gnn = gnn.drop_duplicates()

    gnn.to_csv("GNN_drones" + str(chromosome) + ".csv", index=None)

for chromosome in range(1, 17):
    chrGnn = pd.read_csv("GNN_drones" + str(chromosome) + ".csv", index_col = None)
    if chromosome == 1:
        combinedGnn = chrGnn.iloc[:, 0:4]
    else:
        combinedGnn = combinedGnn + chrGnn.iloc[:, 0:4] * chrPer[chromosome - 1]

combinedGnn['Name'] = [name for (id, name) in sample_ids_names_H]
combinedGnn.to_csv("GNN_hybridDronesCombined.csv", index=None)