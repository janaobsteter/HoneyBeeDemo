import pandas as pd
import scipy
import seaborn as sns
import tskit
import json
import numpy as np
from numpy import genfromtxt
import matplotlib.pyplot as plt
import matplotlib.colors as mc

chromosome = 1

drone_ts = tskit.load("Chr" + str(chromosome) + ".trees")

samples_listed_by_all= [ind.nodes for ind in drone_ts.individuals()]
A = genfromtxt('GNNIndsTable.csv', delimiter=',', skip_header=1)


cols = {drone_ts.node(u).individual: A[:, u] for u in drone_ts.samples()}
cols["population"] = [json.loads(drone_ts.population(drone_ts.node(u).population).metadata)["lineage"] for u in drone_ts.samples()]
cols["individual"] = [drone_ts.node(u).individual for u in drone_ts.samples()]
df = pd.DataFrame(cols)
df.head()

# Plot the data
figsize = (10, 10)
dfg = df.groupby("individual").mean()
# Zscore normalise
for col in list(dfg):
    dfg[col] = scipy.stats.zscore(dfg[col])
row_linkage = scipy.cluster.hierarchy.linkage(dfg, method="average")
order = scipy.cluster.hierarchy.leaves_list(row_linkage)
x_pop = dfg.index.values[order]

cg = sns.clustermap(
    dfg[x_pop], row_linkage=row_linkage, col_cluster=False, # dfg[x_pop] dfgOrder[list(dfgOrder.individual)]
    figsize=figsize, rasterized=True, row_colors = pd.Series(get_ind_colours(df)))


cg.ax_col_dendrogram.legend(ncol=5)
cg.ax_heatmap.set_yticks([])
cg.ax_heatmap.set_xticks([])

for region, col in get_lineage_colours().items():
    cg.ax_col_dendrogram.bar(0, 0, color=col, label=region, linewidth=0)
cg.ax_col_dendrogram.legend(bbox_to_anchor=(50, 0.8))
cg.ax_heatmap.legend(bbox_to_anchor=(0.02, 0.0))
plt.show()

# Try changing the palette
def NonLinCdict(steps, hexcol_array):
    cdict = {'red': (), 'green': (), 'blue': ()}
    for s, hexcol in zip(steps, hexcol_array):
        rgb =matplotlib.colors.hex2color(hexcol)
        cdict['red'] = cdict['red'] + ((s, rgb[0], rgb[0]),)
        cdict['green'] = cdict['green'] + ((s, rgb[1], rgb[1]),)
        cdict['blue'] = cdict['blue'] + ((s, rgb[2], rgb[2]),)
    return cdict

hc = ['#e5e5ff', '#acacdf', '#7272bf', '#39399f', '#000080']
th = [0.000000e+00, 4.868528e-07, 1.648622e-03, 3.296758e-03, 4.944893e-03, 6.593028e-03, 8.241164e-03, 9.889299e-03,
     1.153743e-02, 1.318557e-02, 1.483370e-02, 1.648184e-02, 1.812998e-02, 1.977811e-02, 2.142625e-02, 2.307438e-02,
      2.472252e-02, 2.637065e-02, 2.801879e-02, 1.000000e+00]

cdict = NonLinCdict(th, hc)
cm = mc.LinearSegmentedColormap('test', cdict)


cg = sns.clustermap(
    dfg[x_pop], row_linkage=row_linkage, col_cluster=False, # dfg[x_pop] dfgOrder[list(dfgOrder.individual)]
    figsize=figsize, rasterized=True, row_colors = pd.Series(get_ind_colours(df)))
plt.show()

# COlour by regions
def get_lineage_colours():
    return {
        "A": sns.color_palette("Greens", 2)[1],
        "C": sns.color_palette("Blues", 1)[0],
        "H": sns.color_palette("Greys", 3)[0],
        "M": sns.color_palette("Reds", 2)[1],
        "O": sns.color_palette("Purples", 2)[1],
    }



def get_ind_colours(dataframe):
    # TODO add option to give shades for the different pops.
    lineage_colours = get_lineage_colours()
    ind_colour_map = {}
    for lineage, ind in zip(dataframe.population, dataframe.individual):
        ind_colour_map[ind] = lineage_colours[lineage]
    return ind_colour_map


coloursInd = get_ind_colours()
coloursLineage = get_lineage_colours()
lineage_order = ['A', 'M', 'C', 'O', 'H']

full_df = df
A = np.zeros((len(lineage_order), len(full_df)))



# Try doing just the heatmpa
dfOrder = df.sort_values(by = "population")
dfg = dfOrder.groupby("individual")
dfg = dfg[list(dfOrder.individual)]

Ao = [x[list(dfOrder.individual)] for x in A]
Aoo = [Ao[x] for x in list(dfOrder.individual)]

hm = sns.heatmap(A, vmin=0, vmax=1)