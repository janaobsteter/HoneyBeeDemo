import matplotlib.pyplot as plt
import tskit
import json
import numpy as np
import pandas as pd
import scipy
import seaborn as sns
from scipy import stats
from scipy import cluster


def lineages_colours():
    return {
        "A": sns.color_palette("Oranges", 2)[1],
        "C": sns.color_palette("Purples", 1)[0],
        "H": sns.color_palette("Greys", 3)[0],
        "M": sns.color_palette("Blues", 2)[1],
        "O": sns.color_palette("Greens", 2)[1],
    }

def subspecies_colours():
    return {
        "capensis": "#f7e225",
        "scutellata": "#feb72d",
        "unicolor": "#f07f4f",
        "ligustica": "#310597",
        "carnica": "#ad2793",
        "mellifera": "#228c8d",
        "caucasica": "#a8db34",
        "hybrid": sns.color_palette("Greys", 3)[0]
    }



def get_subspecies_colours(dataframe):
    # TODO add option to give shades for the different pops.
    lineage_colours = lineages_colours()
    subspecies_colour_map = {}
    for lineage, subspecie in zip(dataframe.population, dataframe.subspecie):
        subspecies_colour_map[subspecie] = lineage_colours[lineage]
    return subspecies_colour_map


def get_ind_colours(dataframe):
    # TODO add option to give shades for the different pops.
    lineage_colour = lineages_colours()
    ind_colour_map = {}
    for lineage, ind in zip(dataframe.population, dataframe.individual):
        ind_colour_map[ind] = lineage_colour[lineage]
    return ind_colour_map

def get_ind_colours_subspecies(dataframe):
    # TODO add option to give shades for the different pops.
    subspecie_colours = subspecies_colours()
    ind_colour_map = {}
    for subspecie, ind in zip(dataframe.subspecie, dataframe.individual):
        ind_colour_map[ind] = subspecie_colours[subspecie]
    return ind_colour_map

meta = pd.read_csv("Meta.csv")
countries = list(set(meta.Country))
countries = [x for x in countries if x not in ["USA", ""]]
countries.append("FRAReud")
colourPallete = sns.color_palette("Paired", len(countries))
colourPallete[14] = (0,0,0)
country_colours = dict(zip(countries, colourPallete))

def get_ind_colours_countries(dataframe, country_colours):
    # TODO add option to give shades for the different pops.
    ind_colour_map = {}
    for country, ind in zip(dataframe.country, dataframe.individual):
        ind_colour_map[ind] = country_colours[country]
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


samples_listed_by_all = [ind.nodes for ind in drone_ts.individuals()]

sample_ids_names = [
    (n.id, json.loads(drone_ts.individual(n.individual).metadata)["name"])
    for n in sample_nodes
]

sample_id_country = [
    (id, list(meta.Country[meta.ID == name])[0]) if name in list(meta.ID) else (id, "FRAReud") for id, name in sample_ids_names
]

c = list(set([y for x,y in sample_id_country]))
samples_by_country = [
    np.array([id for (id, name) in sample_id_country if name == country])
    for country in c]


for chromosome in [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]:
    print(chromosome)
    drone_ts = tskit.load("Trees/Chr" + str(chromosome) + "_processed.trees")

    #INDIVIDUAL
    gnn = drone_ts.genealogical_nearest_neighbours(
        drone_ts.samples(), samples_listed_by_all
    )


    cols = {drone_ts.node(u).individual: gnn[:, u] for u in drone_ts.samples()}
    cols["population"] = [json.loads(drone_ts.population(drone_ts.node(u).population).metadata)["lineage"] for u in drone_ts.samples()]
    cols["subspecie"] = [json.loads(drone_ts.population(drone_ts.node(u).population).metadata)["subspecies"] for u in drone_ts.samples()]
    cols["individual"] = [drone_ts.node(u).individual for u in drone_ts.samples()]
    df = pd.DataFrame(cols)
    df.to_csv("Statistics/GnnInd_Chr" + str(chromosome) + ".csv")

    #LINEAGE
    # gnn = drone_ts.genealogical_nearest_neighbours(
    #     drone_ts.samples(), samples_listed_by_population
    # )
    #
    #
    # cols = {json.loads(pop.metadata)['lineage']: gnn[:, pop.id] for pop in drone_ts.populations()}
    # cols["population"] = [json.loads(drone_ts.population(drone_ts.node(u).population).metadata)["lineage"] for u in drone_ts.samples()]
    # cols["individual"] = [drone_ts.node(u).individual for u in drone_ts.samples()]
    # df = pd.DataFrame(cols)
    # nonHybridSamples = [list(samples_listed_by_population[x]) for x in [0,1,2,3,4,6,7]]
    # nonHybridSamples = [x for x in chain.from_iterable(nonHybridSamples)]
    # nonHybridSamples.sort()
    # nonHybridSamples = np.array(nonHybridSamples)
    #
    # HybridSamples = [list(samples_listed_by_population[x]) for x in [5]][0]
    #
    #
    # gnn = drone_ts.genealogical_nearest_neighbours(
    #     drone_ts.samples(), samples_listed_by_population
    # )

    # gnnNonH = drone_ts.genealogical_nearest_neighbours(
    #     nonHybridSamples, samples_listed_by_population_nonH
    # )
    #
    # gnnH = drone_ts.genealogical_nearest_neighbours(
    #     HybridSamples, samples_listed_by_all
    # )

    #COUNTRY
    gnn = drone_ts.genealogical_nearest_neighbours(
        drone_ts.samples(), samples_by_country
    )

    cols = {country: gnn[:, j] for j, country in enumerate(c)}
    cols["population"] = [json.loads(drone_ts.population(drone_ts.node(u).population).metadata)["lineage"] for u in drone_ts.samples()]
    cols["individual"] = [drone_ts.node(u).individual for u in drone_ts.samples()]
    cols["country"] = [country for (id, country) in sample_id_country]
    df = pd.DataFrame(cols)

    df.to_csv("GnnInd_Chr" + str(chromosome) + "_Country.csv")
    chrTree = "chr" + str(chromosome)
    exec(chrTree + " = df")

############################################################################
############################################################################
#Ind
############################################################################
############################################################################
for chromosome in [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]:
    chrGnn = pd.read_csv('Statistics/GnnInd_Chr' + str(chromosome) + '.csv', index_col = None)
    if chromosome == 1:
        combinedGnn = chrGnn.iloc[:, 1:692]
    else:
        combinedGnn = combinedGnn + chrGnn.iloc[:, 1:692] * chrPer[chromosome - 1]
combinedGnn.to_csv("WeightedGnn_Ind.csv", index = False)


combinedGnn = pd.read_csv("WeightedGnn_Ind.csv")
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
lineages = lineages_colours().keys()
lineagesCols = lineages_colours()
handles = [Patch(facecolor=lineagesCols[name]) for name in lineages]
plt.legend(handles, lineagesCols, title='Subspecies',
           bbox_to_anchor=(1, 1), bbox_transform=plt.gcf().transFigure, loc='upper right',
           prop={'size': 20}, ncol=3, title_fontsize=20)
subspecies = subspecies_colours().keys()
subspeciesCols = subspecies_colours()
handles = [Patch(facecolor=subspeciesCols[name]) for name in subspecies]
plt.legend(handles, subspeciesCols, title='Subspecies',
           bbox_to_anchor=(1, 1), bbox_transform=plt.gcf().transFigure, loc='upper right',
           prop={'size': 20}, ncol=3, title_fontsize=20)

plt.show()
cg.figure.savefig("GnnInd_Logit_Lineages.png")
cg.figure.savefig("GnnInd_Logit_SubSpecies.png")


############################################################################
# Only the non-hybrid samples
nonHybridSamples = [list(samples_listed_by_population[x]) for x in [0,1,2,3,4,6,7]]
nonHybridSamples = [item for sublist in nonHybridSamples for item in sublist]

hybridSamples = [list(samples_listed_by_population[x]) for x in [5]]
hybridSamples = [item for sublist in hybridSamples for item in sublist]

combinedGnnGH = combinedGnnG.loc[nonHybridSamples][nonHybridSamples]
row_linkage = scipy.cluster.hierarchy.linkage(combinedGnnGH, method="average")
order = scipy.cluster.hierarchy.leaves_list(row_linkage)
x_pop = combinedGnnGH.index.values[order]
#combinedGnnGH.columns = [int(x) for x in combinedGnnGH.columns]

#Subspecie
cg = sns.clustermap(
    combinedGnnGH[x_pop], row_linkage=row_linkage, col_cluster=False, # dfg[x_pop] dfgOrder[list(dfgOrder.individual)]
    figsize=figsize, rasterized=True,
    row_colors = pd.Series(get_ind_colours_subspecies(combinedGnn[combinedGnn['population']!= "hybrid"])))

cg.ax_heatmap.set_yticks([])
cg.ax_heatmap.set_xticks([])
#cg.cax.set_visible(False)

from matplotlib.patches import Patch

lineages = ['A', 'M', 'C', 'O', 'H']
subspecies = subspecies_colours().keys()
subspeciesCols = subspecies_colours()
handles = [Patch(facecolor=subspeciesCols[name]) for name in subspecies]
plt.legend(handles, subspeciesCols, title='Subspecies',
           bbox_to_anchor=(1, 1), bbox_transform=plt.gcf().transFigure, loc='upper right',
           prop={'size': 20}, ncol=3, title_fontsize=20)

plt.show()

cg.figure.savefig("/home/x/JANA/HoneybeeGeno/tsinfer/GnnInd_Logit_SubSpecies_Hybrid.png")

#country
cg = sns.clustermap(
    combinedGnnGH[x_pop], row_linkage=row_linkage, col_cluster=False, # dfg[x_pop] dfgOrder[list(dfgOrder.individual)]
    figsize=figsize, rasterized=True,
    row_colors = pd.Series(get_ind_colours_countries(combinedGnn[combinedGnn['population']!= "hybrid"], country_colours)))
cg.ax_heatmap.set_yticks([])
cg.ax_heatmap.set_xticks([])

from matplotlib.patches import Patch

countries = [x for x in country_colours.keys() if x in list(combinedGnn.country[combinedGnn['population']== "H"])]
handles = [Patch(facecolor=country_colours[name]) for name in countries]
plt.legend(handles, country_colours, title='Country',
           bbox_to_anchor=(1, 1), bbox_transform=plt.gcf().transFigure, loc='upper right',
           prop={'size': 20}, ncol=4)

plt.show()

cg.figure.savefig("GnnInd_Logit_Countries_Hybrid.png")
############################################################################
# Only the hybrid samples
combinedGnnGnH = combinedGnnG.loc[hybridSamples][hybridSamples]
row_linkage = scipy.cluster.hierarchy.linkage(combinedGnnGnH, method="average")
order = scipy.cluster.hierarchy.leaves_list(row_linkage)
x_pop = combinedGnnGnH.index.values[order]
#combinedGnnGH.columns = [int(x) for x in combinedGnnGH.columns]

#Subspecie
cg = sns.clustermap(
    combinedGnnGnH[x_pop], row_linkage=row_linkage, col_cluster=False, # dfg[x_pop] dfgOrder[list(dfgOrder.individual)]
    figsize=figsize, rasterized=True,
    row_colors = pd.Series(get_ind_colours_countries(combinedGnn[combinedGnn['population']== "H"], country_colours)))

cg.ax_heatmap.set_yticks([])
cg.ax_heatmap.set_xticks([])

from matplotlib.patches import Patch

subspecies = [x for x in country_colours.keys() if x in list(combinedGnn.subspecie[combinedGnn['population']== "H"])]
handles = [Patch(facecolor=country_colours[name]) for name in countries]
plt.legend(handles, country_colours, title='Country',
           bbox_to_anchor=(1, 1), bbox_transform=plt.gcf().transFigure, loc='upper right',
           prop={'size': 20}, ncol=4)

plt.show()


##### Only the FRAReud hybrids###
FRAReud_samples = [x for x, y in sample_id_country if y == "FRAReud"]
FRAReud_nonH_samples = nonHybridSamples + FRAReud_samples

combinedGnnGnHReu = combinedGnnG.loc[FRAReud_nonH_samples][FRAReud_nonH_samples]
row_linkage = scipy.cluster.hierarchy.linkage(combinedGnnGnHReu, method="average")
order = scipy.cluster.hierarchy.leaves_list(row_linkage)
x_pop = combinedGnnGnHReu.index.values[order]

#Subspecie
cg = sns.clustermap(
    combinedGnnGnHReu[x_pop], row_linkage=row_linkage, col_cluster=False, # dfg[x_pop] dfgOrder[list(dfgOrder.individual)]
    figsize=figsize, rasterized=True,
    row_colors = pd.Series(get_ind_colours_subspecies(combinedGnn[combinedGnn['individual'].isin(FRAReud_nonH_samples)])))

#cg.ax_heatmap.set_yticks([])
#cg.ax_heatmap.set_xticks([])

from matplotlib.patches import Patch

subspecies = subspecies_colours().keys()
subspeciesCols = subspecies_colours()
handles = [Patch(facecolor=subspeciesCols[name]) for name in subspecies]
plt.legend(handles, subspeciesCols, title='Subspecies',
           bbox_to_anchor=(1, 1), bbox_transform=plt.gcf().transFigure, loc='upper right',
           prop={'size': 20}, ncol=3, title_fontsize=20, borderaxespad=0)

plt.show()

cg.figure.savefig("GnnInd_NonHybrid_FRAReud.png")

##### Only the ROD hybrids###
ROD_samples = [x for x, y in sample_id_country if y == "ROD"]
ROD_nonH_samples = nonHybridSamples + ROD_samples

combinedGnnGnHROD = combinedGnnG.loc[ROD_nonH_samples][ROD_nonH_samples]
row_linkage = scipy.cluster.hierarchy.linkage(combinedGnnGnHROD, method="average")
order = scipy.cluster.hierarchy.leaves_list(row_linkage)
x_pop = combinedGnnGnHROD.index.values[order]

#Subspecie
cg = sns.clustermap(
    combinedGnnGnHROD[x_pop], row_linkage=row_linkage, col_cluster=False, # dfg[x_pop] dfgOrder[list(dfgOrder.individual)]
    figsize=figsize, rasterized=True,
    row_colors = pd.Series(get_ind_colours_subspecies(combinedGnn[combinedGnn['individual'].isin(ROD_nonH_samples)])))

cg.ax_heatmap.set_yticks([])
cg.ax_heatmap.set_xticks([])

from matplotlib.patches import Patch

subspecies = subspecies_colours().keys()
subspeciesCols = subspecies_colours()
handles = [Patch(facecolor=subspeciesCols[name]) for name in subspecies]
plt.legend(handles, subspeciesCols, title='Subspecies',
           bbox_to_anchor=(1, 1), bbox_transform=plt.gcf().transFigure, loc='upper right',
           prop={'size': 20}, ncol=3, title_fontsize=20)

plt.show()



##### Only the CAN hybrids###
CAN_samples = [x for x, y in sample_id_country if y == "CAN"]
CAN_nonH_samples = nonHybridSamples + CAN_samples

combinedGnnGnHCAN = combinedGnnG.loc[CAN_nonH_samples][CAN_nonH_samples]
row_linkage = scipy.cluster.hierarchy.linkage(combinedGnnGnHCAN, method="average")
order = scipy.cluster.hierarchy.leaves_list(row_linkage)
x_pop = combinedGnnGnHCAN.index.values[order]

#Subspecie
cg = sns.clustermap(
    combinedGnnGnHCAN[x_pop], row_linkage=row_linkage, col_cluster=False, # dfg[x_pop] dfgOrder[list(dfgOrder.individual)]
    figsize=figsize, rasterized=True,
    row_colors = pd.Series(get_ind_colours_subspecies(combinedGnn[combinedGnn['individual'].isin(CAN_nonH_samples)])))

cg.ax_heatmap.set_yticks([])
cg.ax_heatmap.set_xticks([])

from matplotlib.patches import Patch

subspecies = subspecies_colours().keys()
subspeciesCols = subspecies_colours()
handles = [Patch(facecolor=subspeciesCols[name]) for name in subspecies]
plt.legend(handles, subspeciesCols, title='Subspecies',
           bbox_to_anchor=(1, 1), bbox_transform=plt.gcf().transFigure, loc='upper right',
           prop={'size': 20}, ncol=3, title_fontsize=20)

plt.show()
############################################################################
############################################################################
#Country
############################################################################
############################################################################
for chromosome in [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]:
    chrGnn = pd.read_csv('GnnInd_Chr' + str(chromosome) + '_Country.csv', index_col = None)
    if chromosome == 1:
        combinedGnn = chrGnn.iloc[:, 1:16]
    else:
        combinedGnn = combinedGnn + chrGnn.iloc[:, 1:16] * chrPer[chromosome - 1]
combinedGnn.to_csv("WeightedGnn_Ind_Country.csv", index = False)


combinedGnn = pd.read_csv("WeightedGnn_Ind_Country.csv")
combinedGnn["country"] = [country for (id, country) in sample_id_country]
# Plot the data
figsize = (10, 10)
combinedGnnG = combinedGnn.groupby("country").mean()
# combinedGnn['subspecie'] = sample_subspecie
# combinedGnn['individual'] = chrGnn.individual
# combinedGnn['population'] = chrGnn.population


# Zscore normalise
for col in list(combinedGnnG):
    combinedGnnG[col] = [np.log(x+ 0.0001) for x in list(combinedGnnG[col])]
row_linkage = scipy.cluster.hierarchy.linkage(combinedGnnG, method="average")
order = scipy.cluster.hierarchy.leaves_list(row_linkage)
x_pop = combinedGnnG.index.values[order]

#Transform
cg = sns.clustermap(
    combinedGnnG[x_pop], row_linkage=row_linkage, col_cluster=False, # dfg[x_pop] dfgOrder[list(dfgOrder.individual)]
    figsize=figsize, rasterized=True, row_colors = pd.Series(get_ind_colours_subspecies(combinedGnn)))
cg.ax_heatmap.set_yticks([])
cg.ax_heatmap.set_xticks([])
#cg.cax.set_visible(False)

from matplotlib.patches import Patch

lineages = ['A', 'M', 'C', 'O', 'H']
subspecies = subspecies_colours().keys()
subspeciesCols = subspecies_colours()
handles = [Patch(facecolor=subspeciesCols[name]) for name in subspecies]
plt.legend(handles, subspeciesCols, title='Subspecies',
           bbox_to_anchor=(1, 1), bbox_transform=plt.gcf().transFigure, loc='upper right',
           prop={'size': 20}, ncol=3)

plt.show()


# Only the non-hybrid samples
nonHybridSamples = [list(samples_listed_by_population[x]) for x in [0,1,2,3,4,6,7]]
nonHybridSamples = [item for sublist in nonHybridSamples for item in sublist]

hybridSamples = [list(samples_listed_by_population[x]) for x in [5]]
hybridSamples = [item for sublist in hybridSamples for item in sublist]


combinedGnnH = combinedGnn.loc[nonHybridSamples]
combinedGnnGH = combinedGnnH.groupby("country").mean()
for col in list(combinedGnnGH):
    combinedGnnGH[col] = [np.log(x+ 0.0001) for x in list(combinedGnnGH[col])]
row_linkage = scipy.cluster.hierarchy.linkage(combinedGnnGH, method="average")
order = scipy.cluster.hierarchy.leaves_list(row_linkage)
x_pop = combinedGnnGH.index.values[order]
#combinedGnnGH.columns = [int(x) for x in combinedGnnGH.columns]

#Subspecie
cg = sns.clustermap(
    combinedGnnGH, row_linkage=row_linkage, col_cluster=False, # dfg[x_pop] dfgOrder[list(dfgOrder.individual)]
    figsize=figsize, rasterized=True,
    row_colors = pd.Series(get_ind_colours_subspecies(combinedGnn[combinedGnn['population']!= "hybrid"])))

cg.ax_heatmap.set_yticks([])
cg.ax_heatmap.set_xticks([])
#cg.cax.set_visible(False)

from matplotlib.patches import Patch

lineages = ['A', 'M', 'C', 'O', 'H']
subspecies = subspecies_colours().keys()
subspeciesCols = subspecies_colours()
handles = [Patch(facecolor=subspeciesCols[name]) for name in subspecies]
plt.legend(handles, subspeciesCols, title='Subspecies',
           bbox_to_anchor=(1, 1), bbox_transform=plt.gcf().transFigure, loc='upper right',
           prop={'size': 20}, ncol=3)

plt.show()

#country
cg = sns.clustermap(
    combinedGnnGH, row_linkage=row_linkage, col_cluster=False, # dfg[x_pop] dfgOrder[list(dfgOrder.individual)]
    figsize=figsize, rasterized=True,
    row_colors = pd.Series(get_ind_colours_countries(combinedGnn[combinedGnn['population']!= "hybrid"], country_colours)))
cg.ax_heatmap.set_yticks([])
cg.ax_heatmap.set_xticks([])
#cg.cax.set_visible(False)

from matplotlib.patches import Patch

countries = country_colours.keys()
handles = [Patch(facecolor=country_colours[name]) for name in countries]
plt.legend(handles, country_colours, title='Country',
           bbox_to_anchor=(1, 1), bbox_transform=plt.gcf().transFigure, loc='upper right',
           prop={'size': 20}, ncol=5)

plt.show()



############################################################################
############################################################################
#Subspecie
############################################################################
############################################################################
for chromosome in [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]:
    chrGnn = pd.read_csv('Statistics/GnnInd_Chr' + str(chromosome) + '_Subspecies.csv', index_col = None)
    if chromosome == 1:
        combinedGnn = chrGnn.iloc[:, 1:9]
    else:
        combinedGnn = combinedGnn + chrGnn.iloc[:, 1:9] * chrPer[chromosome - 1]

combinedGnn['subspecie'] = sample_subspecie
# combinedGnn['subspecie'] = [sample_subspecie[x] for x in nonHybridSamples]
# combinedGnn['subspecie'] = [sample_subspecie[x] for x in HybridSamples]
combinedGnn['individual'] = chrGnn.individual
combinedGnn['population'] = chrGnn.population



# Plot the data
figsize = (10, 10)
combinedGnnG = combinedGnn.groupby("subspecie").mean()
# Zscore normalise
for col in list(combinedGnnG):
    combinedGnnG[col] = scipy.stats.zscore(combinedGnnG[col])
row_linkage = scipy.cluster.hierarchy.linkage(combinedGnnG, method="average")
order = scipy.cluster.hierarchy.leaves_list(row_linkage)
x_pop = combinedGnnG.index.values[order]

cg = sns.clustermap(
    combinedGnnG[x_pop], row_linkage=row_linkage, col_cluster=False, # dfg[x_pop] dfgOrder[list(dfgOrder.individual)]
    figsize=figsize, rasterized=True, row_colors = pd.Series(get_subspecies_colours(combinedGnn)))
plt.show()

############################################################################
############################################################################
# Lineage
############################################################################
############################################################################
for chromosome in [1,2,3,4,5,6,8,9,10,11,12,13,14,15,16]:
    chrGnn = pd.read_csv('GnnInd_Chr' + str(chromosome) + '_Subspecies.csv', index_col = None)
    if chromosome == 1:
        combinedGnn = chrGnn.iloc[:, 1:9]
    else:
        combinedGnn = combinedGnn + chrGnn.iloc[:, 1:9] * chrPer[chromosome - 1]

#combinedGnn['population'] = chrGnn.population
combinedGnn['subspecie'] = sample_subspecie
# Plot the data
figsize = (10, 10)
combinedGnnG = combinedGnn.groupby("subspecie").mean()
# Zscore normalise
for col in list(combinedGnnG):
    combinedGnnG[col] = scipy.stats.zscore(combinedGnnG[col])
row_linkage = scipy.cluster.hierarchy.linkage(combinedGnnG, method="average")
order = scipy.cluster.hierarchy.leaves_list(row_linkage)
x_pop = combinedGnnG.index.values[order]

cg = sns.clustermap(
    combinedGnnG[x_pop], row_linkage=row_linkage, col_cluster=False, # dfg[x_pop] dfgOrder[list(dfgOrder.individual)]
    figsize=figsize, rasterized=True, row_colors = pd.Series(get_lineage_colours()))

plt.show()

# Subspecie
for chromosome in [1,2,3,4,5,6,8,9,10,11,12,13,14,15,16]:
    chrGnn = pd.read_csv('Statistics/GnnInd_Chr' + str(chromosome) + '_Lineage.csv', index_col = None)
    if chromosome == 1:
        combinedGnn = chrGnn.iloc[:, 1:6]
    else:
        combinedGnn = combinedGnn + chrGnn.iloc[:, 1:6] * chrPer[chromosome - 1]

combinedGnn['population'] = chrGnn.population
# Plot the data
figsize = (10, 10)
combinedGnnG = combinedGnn.groupby("population").mean()
# Zscore normalise
for col in list(combinedGnnG):
    combinedGnnG[col] = scipy.stats.zscore(combinedGnnG[col])
row_linkage = scipy.cluster.hierarchy.linkage(combinedGnnG, method="average")
order = scipy.cluster.hierarchy.leaves_list(row_linkage)
x_pop = combinedGnnG.index.values[order]

cg = sns.clustermap(
    combinedGnnG[x_pop], row_linkage=row_linkage, col_cluster=False, # dfg[x_pop] dfgOrder[list(dfgOrder.individual)]
    figsize=figsize, rasterized=True, row_colors = pd.Series(get_subspecies_colours()))


# By Ind
gnn = drone_ts.genealogical_nearest_neighbours(
    drone_ts.samples(), samples_listed_by_all
)



columns = [json.loads(ind.metadata)["name"] for ind in drone_ts.individuals()]
gnn_table = pd.DataFrame(
    data=gnn,
    index=[
        pd.Index(sample_ids, name="Sample node"),
        pd.Index(sample_names, name="Drones"),
        pd.Index(sample_pops, name="Lineage"),
    ],
    columns=columns,
)


cols = {json.loads(pop.metadata)["lineage"]: gnn[:, pop.id] for pop in drone_ts.populations()}
cols["population"] = [json.loads(drone_ts.population(drone_ts.node(u).population).metadata)["lineage"] for u in drone_ts.samples()]
cols["individual"] = [drone_ts.node(u).individual for u in drone_ts.samples()]
df = pd.DataFrame(cols)