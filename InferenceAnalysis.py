import tskit
import tsinfer
import json
import pandas as pd
import numpy as np
import seaborn as sns
import scipy

chromosome = 3

drone_ts = tskit.load("Chr" + str(chromosome) + ".trees")

print("Draw a tree")
# Set the colours of the subspecies for plotting
colours = {"hybrid": "grey", "carnica": "yellow", "ligustica": "orange", "mellifera": "red", "caucasica": "blue",
           "scutellata": "black", "unicolor": "blue", "capensis": "purple", }
# Set the colours for the individuals
colours_for_node = {}

# If you want to plot by lineages
# Set the colours of the lineages for plotting
coloursL = {"A": "black", "C": "yellow", "M": "red", "O": "blue", "H": "grey"}
colours_for_nodeL = {}
for n in drone_ts.samples():
    population_data = drone_ts.population(drone_ts.node(n).population)
    colours_for_node[n] = colours[json.loads(population_data.metadata)["subspecies"]]
    colours_for_nodeL[n] = coloursL[json.loads(population_data.metadata)["lineage"]]

# Get the sample names for the sample nodes
individual_for_node = {}
for n in drone_ts.samples():
    individual_data = drone_ts.individual(drone_ts.node(n).individual)
    individual_for_node[n] = json.loads(individual_data.metadata)["name"]

# Obtain node info - but just for samples!
# There is also intermediate nodes!!! (total #nodes > #nodes for samples)
sample_nodes = [drone_ts.node(n) for n in drone_ts.samples()]

# Get samples ids
sample_ids = [n.id for n in sample_nodes]

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

samples_listed_by_population = [
    drone_ts.samples(population=pop_id)
    for pop_id in range(drone_ts.num_populations)
]


samples_listed_by_all= [ind.nodes for ind in drone_ts.individuals()]
# This one takes a long time since there is many nodes!!!
# node_names = [
#     json.loads(drone_ts.individual(n.individual).metadata)["name"]
#     if n in sample_nodes else None for n in drone_ts.nodes()
# ]

# Create a dictionary holding the sample node ids within each population (in this case sample_pops are lineages)
pops = dict()
for (node, pop) in zip(sample_ids, sample_pops):
    if not pop in pops.keys():
        pops[pop] = [node]
    else:
        pops[pop].append(node)

print("Do the statistics.")
# Obtain Fst and write a table
Fst = drone_ts.Fst([pops["A"], pops["C"], pops["M"], pops["O"]],
                   indexes=[(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)])
# print(Fst)
fstDF = pd.DataFrame(Fst, columns=["Fst"])
fstDF.loc[:, "Group1"] = ["A", "A", "A", "C", "C", "M"]
fstDF.loc[:, "Group2"] = ["C", "M", "O", "M", "O", "O"]
fstDF.to_csv("FST" + str(chromosome) + ".csv", index=None)

# Obtain Tajima's D and write a table
tajimaD = drone_ts.Tajimas_D([pops["A"], pops["C"], pops["M"], pops["O"], pops["H"]])
# print(tajimaD)
tajimaDF = pd.DataFrame(tajimaD, columns=["TajimasD"])
tajimaDF.loc[:, "Group"] = ["A", "C", "M", "O", "H"]
tajimaDF.to_csv("TajimasD" + str(chromosome) + ".csv", index=None)

# Divergence
divergence = drone_ts.divergence([pops["A"], pops["C"], pops["M"], pops["O"]],
                                 indexes=[(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)])
# print(divergence)
divergenceDF = pd.DataFrame(divergence, columns=["Divergence"])
divergenceDF.loc[:, "Group1"] = ["A", "A", "A", "C", "C", "M"]
divergenceDF.loc[:, "Group2"] = ["C", "M", "O", "M", "O", "O"]
divergenceDF.to_csv("Divergence" + str(chromosome) + ".csv", index=None)

# Diversity and write a table
diversity = drone_ts.diversity([pops["A"], pops["C"], pops["M"], pops["O"], pops["H"]])
diversityDF = pd.DataFrame(diversity, columns=["TajimasD"])
diversityDF.loc[:, "Group"] = ["A", "C", "M", "O", "H"]
diversityDF.to_csv("Diversity" + str(chromosome) + ".csv", index=None)

# Allele frequency spectrum
np.savetxt("AFS_Alineage" + str(chromosome) + ".csv", drone_ts.allele_frequency_spectrum([pops["A"]]))
np.savetxt("AFS_Mlineage" + str(chromosome) + ".csv", drone_ts.allele_frequency_spectrum([pops["M"]]))
np.savetxt("AFS_Clineage" + str(chromosome) + ".csv", drone_ts.allele_frequency_spectrum([pops["C"]]))
np.savetxt("AFS_Olineage" + str(chromosome) + ".csv", drone_ts.allele_frequency_spectrum([pops["O"]]))

# f2, f3, f4 statistics
f2 = drone_ts.f2([pops["A"], pops["C"], pops["M"], pops["O"]], indexes=[(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)])
# print(f2)
f2DF = pd.DataFrame(f2, columns=["f2"])
f2DF.loc[:, "Group1"] = ["A", "A", "A", "C", "C", "M"]
f2DF.loc[:, "Group2"] = ["C", "M", "O", "M", "O", "O"]
f2DF.to_csv("f2" + str(chromosome) + ".csv", index=None)

# f4
f4 = drone_ts.f4([pops["A"], pops["C"], pops["M"], pops["O"]])
print("This is f4 statistics: ", str(f4))

# GNN
gnn = pd.DataFrame(drone_ts.genealogical_nearest_neighbours(pops["H"], sample_sets=[pops["A"], pops["C"], pops["M"], pops["O"]]), columns=["A", "C", "M", "O"])

# Add sample ids and names to the data frame
# gnn.loc[:, "Id"] = pops["H"] #The two nodes per individual will be the same since they are homozygous diploids
namesH = [name for (name, ID) in zip(sample_names, sample_ids) if ID in pops['H']]
gnn.loc[:, "Names"] = namesH
gnn = gnn.drop_duplicates()

gnn.to_csv("GNN_drones" + str(chromosome) + ".csv", index=None)

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


gnn_table.to_csv("GNNIndsTable_chr2.csv", index=False, header=True)

def get_tgp_colours():
    # TODO add option to give shades for the different pops.
    region_colours = get_tgp_region_colours()
    pop_colour_map = {}
    for region, pops in tgp_region_pop.items():
        for pop in pops:
            pop_colour_map[pop] = region_colours[region]
    return pop_colour_map

def get_tgp_region_colours():
    return {
        "A": sns.color_palette("Greens", 2)[1],
        "C": sns.color_palette("Blues", 1)[0],
        "H": sns.color_palette("Greys", 3)[0],
        "M": sns.color_palette("Reds", 2)[1],
        "O": sns.color_palette("Purples", 2)[1],
    }


tgp_region_pop = {
    'AMR': ['CLM', 'MXL', 'PUR', 'PEL'],
    'AFR': ['LWK', 'ASW', 'GWD', 'MSL', 'YRI', 'ACB', 'ESN'],
    'EAS': ['CHS', 'KHV', 'JPT', 'CHB', 'CDX'],
    'SAS': ['BEB', 'STU', 'GIH', 'PJL', 'ITU'],
    'EUR': ['FIN', 'GBR', 'IBS', 'CEU', 'TSI']
}


#df = pd.read_csv("data/1kg_gnn.csv")
dfg = gnn_table.groupby("Lineage").mean()
# Zscore normalise
for col in list(dfg):
    dfg[col] = scipy.stats.zscore(dfg[col])

row_linkage = scipy.cluster.hierarchy.linkage(dfg, method="average")
order = scipy.cluster.hierarchy.leaves_list(row_linkage)
x_pop = dfg.index.values[order]

colours = pd.Series(get_tgp_region_colours())
cg = sns.clustermap(
    dfg[x_pop], row_linkage=row_linkage, col_cluster=False,
    row_colors=colours, figsize=figsize, rasterized=True)
cg.ax_heatmap.set_ylabel("")

for region, col in region_colours.items():
    cg.ax_col_dendrogram.bar(0, 0, color=col, label=region, linewidth=0)
return cg

gnn_table = pd.DataFrame(
    data=gnn,
    index=[
        pd.Index(sample_ids, name="Sample node"),
        pd.Index(sample_names, name="Drones"),
        pd.Index(sample_pops, name="Lineage"),
    ],
    columns=[p.id for p in drone_ts.populations()],
)

# Genetic relatedness
grel = drone_ts.genetic_relatedness([pops["A"], pops["C"], pops["M"], pops["O"]],
                                    indexes=[(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)])
# print(grel)
grelDF = pd.DataFrame(grel, columns=["GeneticRel"])
grelDF.loc[:, "Group1"] = ["A", "A", "A", "C", "C", "M"]
grelDF.loc[:, "Group2"] = ["C", "M", "O", "M", "O", "O"]
grelDF.to_csv("GeneticRelatedness" + str(chromosome) + ".csv", index=None)

# Mean descendants of ALL nodes!!! Also intermediate
meanDesc = pd.DataFrame(drone_ts.mean_descendants([pops["A"], pops["C"], pops["M"], pops["O"]]),
                        columns=["A", "C", "M", "O"])
meanDesc.loc[:, "NodeID"] = [n.id for n in drone_ts.nodes()]
meanDesc.loc[:, "NodeName"] = None

# Add names for sample nodes
for (ID, name) in zip(sample_nodes, sample_names):
    meanDesc.at[meanDesc.NodeID == ID.id, "NodeName"] = name

meanDesc.to_csv("MeanDescendants" + str(chromosome) + ".csv", index=None)

# Write the tables with mutations, nodes, sites, populations and individuals
drone_ts.dump_text(mutations=open("Mutations" + str(chromosome) + ".txt", "w"),
                   nodes=open("Nodes" + str(chromosome) + ".txt", "w"),
                   sites=open("Sites" + str(chromosome) + ".txt", "w"),
                   edges=open("Edges" + str(chromosome) + ".txt", "w"),
                   individuals=open("Individuals" + str(chromosome) + ".txt", "w"),
                   populations=open("Populations" + str(chromosome) + ".txt", "w", encoding="utf-8"))
