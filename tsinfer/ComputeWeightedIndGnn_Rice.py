import matplotlib.pyplot as plt
import tskit
import json
import numpy as np
import pandas as pd
import scipy
import seaborn as sns
from scipy import stats
from scipy import cluster
from matplotlib.patches import Patch

def subspecies_colours():
    return {
        "japonice": sns.color_palette("Oranges", 2)[1],
        "indica": sns.color_palette("Purples", 1)[0]
    }

def get_ind_colours(dataframe):
    subspecies_colour = subspecies_colours()
    ind_colour_map = {}
    for pop, ind in zip(dataframe.subspecies, dataframe.individual):
        ind_colour_map[ind] = subspecies_colour[pop]
    return ind_colour_map

# Here you will have to input rice data (I have this from stdpopsim, hence the form - the only thing you need is the lengths)
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

# Here you compute the weights based on the chromosom length
chromosomes = genome['chromosomes']
chrLengths = [value['length'] for value in chromosomes.values()]
chrPer = [x / sum(chrLengths) for x in chrLengths]



# Obtain node info - but just for samples! (here you read in only one, the smallest, to get out the sample ids and names)
rice_ts = tskit.load("Trees/Chr16_processed.trees")

#Here you prepare the list for the GNN
sample_nodes = [rice_ts.node(n) for n in rice_ts.samples()]
samples_listed_by_all = [ind.nodes for ind in rice_ts.individuals()]

# Change to however chromosomes you have
for chromosome in range(1, 17):
    print(chromosome)
    # Here set the paths to the trees
    rice_ts = tskit.load("Trees/Chr" + str(chromosome) + "_processed.trees")

    #INDIVIDUAL
    gnn = rice_ts.genealogical_nearest_neighbours(
        rice_ts.samples(), samples_listed_by_all
    )


    cols = {rice_ts.node(u).individual: gnn[:, u] for u in rice_ts.samples()}
    cols["subspecies"] = [json.loads(rice_ts.subspecies(rice_ts.node(u).subspecies).metadata)["subspecies"] for u in rice_ts.samples()]
    cols["individual"] = [rice_ts.node(u).individual for u in rice_ts.samples()]
    df = pd.DataFrame(cols)
    # Write to wherever you want
    df.to_csv("Statistics/GnnInd_Chr" + str(chromosome) + ".csv")



############################################################################
# This here now reads them back in and combines the results for all the chromosomes (weights them)
############################################################################
for chromosome in range(1, 17):
    chrGnn = pd.read_csv('Statistics/GnnInd_Chr' + str(chromosome) + '.csv', index_col = None)
    if chromosome == 1:
        combinedGnn = chrGnn.iloc[:, 1:692]
    else:
        combinedGnn = combinedGnn + chrGnn.iloc[:, 1:692] * chrPer[chromosome - 1]

combinedGnn.to_csv("WeightedGnn_Ind.csv", index = False)


############################################################################
# Plot the data
############################################################################
combinedGnn = pd.read_csv("WeightedGnn_Ind.csv")
# Plot the data
figsize = (10, 10)
combinedGnnG = combinedGnn.copy() #.groupby("individual").mean()
combinedGnn['individual'] = chrGnn.individual
combinedGnn['subspecies'] = chrGnn.subspecies
# Logit normalise (they used zscore normalization - but that didn't work for me!)
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


# This adds the legend
subspecies = subspecies_colours().keys()
subspeciesCols = subspecies_colours()
handles = [Patch(facecolor = subspeciesCols[name]) for name in subspecies]
plt.legend(handles, subspeciesCols, title='Subspecies',
           bbox_to_anchor=(1, 1), bbox_transform=plt.gcf().transFigure, loc='upper right',
           prop={'size': 20}, ncol=3, title_fontsize=20)
plt.show()
cg.figure.savefig("GnnInd_Logit_Subspecies.png")


