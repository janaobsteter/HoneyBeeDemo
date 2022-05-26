import pandas as pd
import scipy
import seaborn as sns
import tskit
import json
import matplotlib.pyplot as plt
import numpy as np


chromosome = 1

lineage_subspecies = {
    'A': ['scutellata', 'unicolor', 'capensis'],
    'C': ['carnica', 'ligustica'],
    'M': ['mellifera'],
    'O': ['caucasica'],
    'H': ['hybrid']
}

subspecies = ['scutellata', 'unicolor', 'capensis', 'carnica', 'ligustica', 'mellifera', 'caucasica', 'hybrid']

# COlour by regions
def get_lineage_colours():
    return {
        "A": sns.color_palette("Greens", 2)[1],
        "C": sns.color_palette("Blues", 1)[0],
        "H": sns.color_palette("Greys", 3)[0],
        "M": sns.color_palette("Reds", 2)[1],
        "O": sns.color_palette("Purples", 2)[1],
    }



drone_ts = tskit.load("Chr" + str(chromosome) + ".trees")

sample_nodes = [drone_ts.node(n) for n in drone_ts.samples()]
sample_subspecie = [
    json.loads(drone_ts.population(n.population).metadata)["subspecies"]
    for n in sample_nodes
]
samples_listed_by_subspecie = [
    drone_ts.samples(population=pop_id)
    for pop_id in range(drone_ts.num_populations)
]
A = drone_ts.genealogical_nearest_neighbours(
    drone_ts.samples(), samples_listed_by_subspecie, num_threads=1)

cols = {json.loads(pop.metadata)["subspecies"]: A[:, pop.id] for pop in drone_ts.populations()}
cols["lineage"] = [json.loads(drone_ts.population(drone_ts.node(u).population).metadata)["lineage"] for u in drone_ts.samples()]
cols["subspecie"] = [json.loads(drone_ts.population(drone_ts.node(u).population).metadata)["subspecies"] for u in drone_ts.samples()]
cols["individual"] = [drone_ts.node(u).individual for u in drone_ts.samples()]
df = pd.DataFrame(cols)
df.head()


pop_colours = {}
lineage_colours = get_lineage_colours()
for subspecie, lineage in df.groupby(["subspecie", "lineage"]).size().index:
    pop_colours[subspecie] = lineage_colours[lineage]


dfg = df.groupby("subspecie").mean()
# Zscore normalise
for col in list(dfg):
    dfg[col] = scipy.stats.zscore(dfg[col])

row_linkage = scipy.cluster.hierarchy.linkage(dfg, method="average")
order = scipy.cluster.hierarchy.leaves_list(row_linkage)
x_pop = dfg.index.values[order]

colours = pd.Series(pop_colours)
cg = sns.clustermap(
    dfg[x_pop], row_linkage=row_linkage, col_cluster=False,
    row_colors=colours, figsize=(10,10), rasterized=True)

cg.ax_heatmap.set_ylabel("")
cg.ax_col_dendrogram.legend(ncol=8)
for lineage, col in lineage_colours.items():
    cg.ax_col_dendrogram.bar(0, 0, color=col, label=lineage, linewidth=0)

from matplotlib.patches import Patch

lineages = ['A', 'M', 'C', 'O', 'H']
handles = [Patch(facecolor=lineage_colours[name]) for name in lineages]
plt.legend(handles, lineage_colours, title='Lineages',
           bbox_to_anchor=(1, 1), bbox_transform=plt.gcf().transFigure, loc='upper right',
           prop={'size': 20}, ncol=8)
#cg.ax_heatmap.set_yticklabels(cg.ax_heatmap.get_yticklabels(), fontsize=20)
#g.ax_heatmap.set_xticks([])

plt.show()

########################################################3
# Plot by country
########################################################3
meta = pd.read_csv("Meta.csv")
sample_names = [
    json.loads(drone_ts.individual(n.individual).metadata)["name"]
    for n in sample_nodes
]

sample_ids_names = [
    (n.id, json.loads(drone_ts.individual(n.individual).metadata)["name"])
    for n in sample_nodes
]
sampleIDNamesDF = pd.DataFrame({"ID": [id for id, name in sample_ids_names],
                                "Name": [name for id, name in sample_ids_names]})
sampleIDNamesDF.to_csv("~/JANA/HoneybeeGeno/tsinfer/SampleIdNames.csv", index=False)
sample_id_country = [
    (id, list(meta.Country[meta.ID == name])[0]) if name in list(meta.ID) else (id, "FRAReud") for id, name in sample_ids_names
]

countries = list(set(meta.Country)) + [""]
countries = [x for x in countries if x not in ["USA", ""]]
countries.append("FRAReud")
common = set(meta.ID).intersection(set(sample_names))
namesNotCountry = [x for x in sample_names if x not in list(common)]
len(namesNotCountry)
namesNotCountry


samples_by_country = [
    np.array([id for (id, name) in sample_ids_names if name in list(meta.ID[meta.Country == country])])
    for country in countries
]
samples_by_country[len(samples_by_country) -1] = np.array([id for (id, name) in sample_ids_names if name in namesNotCountry])

A = drone_ts.genealogical_nearest_neighbours(
    drone_ts.samples(), samples_by_country, num_threads=1)
A = genfromtxt('GNNIndsTable.csv', delimiter=',', skip_header=1)

#Country to country
cols = {country: A[:, j] for j, country in enumerate(countries)}
# Indicidual to country
cols = {drone_ts.node(u).individual: A[:, u] for u in drone_ts.samples()}
cols["population"] = [json.loads(drone_ts.population(drone_ts.node(u).population).metadata)["lineage"] for u in drone_ts.samples()]
cols["individual"] = [drone_ts.node(u).individual for u in drone_ts.samples()]
cols["country"] = [country for (id, country) in sample_id_country]

df = pd.DataFrame(cols)


dfg = df.groupby("individual").mean()
# Zscore normalise
for col in list(dfg):
    dfg[col] = scipy.stats.zscore(dfg[col])

row_linkage = scipy.cluster.hierarchy.linkage(dfg, method="average")
order = scipy.cluster.hierarchy.leaves_list(row_linkage)
x_pop = dfg.index.values[order]

colourPallete = sns.color_palette("Paired", len(countries))
colourPallete[14] = (0,0,0)
country_colours = dict(zip(countries, colourPallete))

def get_ind_colours_Countries(dataframe, country_colours):
    # TODO add option to give shades for the different pops.
    ind_colour_map = {}
    for country, ind in zip(dataframe.country, dataframe.individual):
        ind_colour_map[ind] = country_colours[country]
    return ind_colour_map


cg = sns.clustermap(
    dfg[x_pop], row_linkage=row_linkage, col_cluster=False,
    row_colors=pd.Series(get_ind_colours_Countries(df, country_colours)), figsize=(8, 8), rasterized=True)
cg.ax_heatmap.set_ylabel("")
cg.ax_heatmap.set_yticks([])
cg.ax_heatmap.set_xticks([])

handles = [Patch(facecolor=country_colours[country]) for country in countries]
plt.legend(handles, country_colours, title='Countries',
           bbox_to_anchor=(1, 1), bbox_transform=plt.gcf().transFigure, loc='upper right',
           prop={'size': 8}, ncol=8)


plt.show()



# Try plotting the heatmap
dfgOrder =