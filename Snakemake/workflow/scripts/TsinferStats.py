import tsinfer
import pandas as pd

treeFile = snakemake.input[0]
fstOutput = snakemake.output[0]
tajimasdOutput = snakemake.output[1]
divergenceOutput = snakemake.output[2]
diversityOutput = snakemake.output[3]
f2Output = snakemake.output[4]
gRelOutput = snakemake.output[5]
meanDescOutput = snakemake.output[6]
mutationsOutput = snakemake.output[7]
nodesOutput = snakemake.output[8]
sitesOutput = snakemake.output[9]
edgesOutput = snakemake.output[10]
individualsOutput = snakemake.output[11]
populationsOutput = snakemake.output[12]

ts = tsinfer.load(treeFile)

# Get the sample names for the sample nodes
individual_for_node = {}
for n in ts.samples():
    individual_data = ts.individual(ts.node(n).individual)
    individual_for_node[n] = json.loads(individual_data.metadata)["name"]

#Obtain node info - but just for samples!
# There is also intermediate nodes!!! (total #nodes > #nodes for samples)
sample_nodes = [ts.node(n) for n in ts.samples()]

#Get samples ids
sample_ids = [n.id for n in sample_nodes]

# Get sample names
sample_names = [
    json.loads(ts.individual(n.individual).metadata)["name"]
    for n in sample_nodes
]
# Get sample population
sample_pops = [
    json.loads(ts.population(n.population).metadata)["lineage"]
    for n in sample_nodes
]

# This one takes a long time since there is many nodes!!!
# node_names = [
#     json.loads(ts.individual(n.individual).metadata)["name"]
#     if n in sample_nodes else None for n in ts.nodes()
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
indeces = [(x,y) for x in range(len(pops)) for y in range(len(pops))]
Fst = ts.Fst([pops.values()], indexes=indeces)
#print(Fst)
fstDF = pd.DataFrame(Fst, columns = ["Fst"])
fstDF.loc[:, "Group1"] = [list(pops.keys())[x] for x in [x for (x,y) in indeces]]
fstDF.loc[:, "Group2"] = [list(pops.keys())[x] for x in [y for (x,y) in indeces]]
fstDF.to_csv(fstOutput, index=None)

# Obtain Tajima's D and write a table
tajimaD = ts.Tajimas_D([pops.values()])
#print(tajimaD)
tajimaDF = pd.DataFrame(tajimaD, columns = ["TajimasD"])
tajimaDF.loc[:, "Group"] = list(pops.keys())
tajimaDF.to_csv(tajimasdOutput, index=None)

# Divergence
divergence = ts.divergence([pops.values()], indexes=indeces)
#print(divergence)
divergenceDF = pd.DataFrame(divergence, columns = ["Divergence"])
divergenceDF.loc[:, "Group1"] = [list(pops.keys())[x] for x in [x for (x,y) in indeces]]
divergenceDF.loc[:, "Group2"] = [list(pops.keys())[x] for x in [y for (x,y) in indeces]]
divergenceDF.to_csv(divergenceOutput, index=None)

# Diversity and write a table
diversity = ts.diversity([pops.values()])
diversityDF = pd.DataFrame(diversity, columns = ["Diversity"])
diversityDF.loc[:, "Group"] = list(pops.keys())
diversityDF.to_csv(diversityOutput, index=None)
#
# # Allele frequency spectrum
# for pop, nodes in list(pop.items()):
#     np.savetxt("AFS_" + pop + "lineage" + str(chromosome) + ".csv", ts.allele_frequency_spectrum([nodes]))


# f2, f3, f4 statistics
f2 = ts.f2([pops.values()], indexes=indeces)
f2DF = pd.DataFrame(f2, columns = ["f2"])
f2DF.loc[:, "Group1"] = [list(pops.keys())[x] for x in [x for (x,y) in indeces]]
f2DF.loc[:, "Group2"] = [list(pops.keys())[x] for x in [y for (x,y) in indeces]]
f2DF.to_csv(f2Output, index=None)


# # GNN
# gnn = pd.DataFrame(ts.genealogical_nearest_neighbours(pops["H"], sample_sets=[pops["A"], pops["C"], pops["M"], pops["O"]]), columns = ["A", "C","M", "O"])
#
# # Add sample ids and names to the data frame
# #gnn.loc[:, "Id"] = pops["H"] #The two nodes per individual will be the same since they are homozygous diploids
# namesH = [name for (name, ID) in zip(sample_names, sample_ids) if ID in pops['H']]
# gnn.loc[:, "Names"] = namesH
# gnn = gnn.drop_duplicates()
#
# gnn.to_csv("GNN_drones" + str(chromosome) + ".csv", index=None)

# Genetic relatedness
grel = ts.genetic_relatedness([pops.values()], indexes=indeces)
#print(grel)
grelDF = pd.DataFrame(grel, columns = ["GeneticRel"])
grelDF.loc[:, "Group1"] = [list(pops.keys())[x] for x in [x for (x,y) in indeces]]
grelDF.loc[:, "Group2"] = [list(pops.keys())[x] for x in [y for (x,y) in indeces]]
grelDF.to_csv(gRelOutput, index=None)

# Mean descendants of ALL nodes!!! Also intermediate
meanDesc = pd.DataFrame(ts.mean_descendants([pops.values()]), columns = list(pops.keys()))
meanDesc.loc[:, "NodeID"] = [n.id for n in ts.nodes()]
meanDesc.loc[:, "NodeName"] = None
# Add names for sample nodes
for (ID, name) in zip(sample_nodes, sample_names):
    meanDesc.at[meanDesc.NodeID == ID.id, "NodeName"] = name
meanDesc.to_csv(meanDescOutput, index=None)


# Write the tables with mutations, nodes, sites, populations and individuals
ts.dump_text(mutations=open(mutationsOutput, "w"),
             nodes = open(nodesOutput, "w"),
             sites = open(sitesOutput, "w"),
             edges = open(edgesOutput, "w"),
             encoding="utf-8")
