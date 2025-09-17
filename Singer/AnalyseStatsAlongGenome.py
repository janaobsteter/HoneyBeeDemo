# %%
import tskit
import pandas as pd
import numpy as np
from lets_plot import *
LetsPlot.setup_html()
from collections import Counter
import sys

chromosome = sys.argv[1]


# %%

chromosome = sys.argv[1]

# %%
ts = tskit.load("/sietch_colab/jobste/Singer/Chr" + str(chromosome) + "/Chr" + str(chromosome) + "_NonHybrid_99.trees")

# %%
meta = pd.read_csv("/home/jobste/Data/TsinferMetaData_Vcf.txt")
meta.head()
allMeta = pd.read_csv("/home/jobste/Data/DataDescriptor_SampleMeta.csv")
allMeta.head()
meta.loc[:, "Country"] = [allMeta.Country[allMeta.Run == x].to_list()[0] for x in meta.ID]

samples = pd.read_csv("/home/jobste/Data/Remapped/NonHybrid_samples.txt", header = None, names = ["ID"])

# %%
# Add subpop as pop

sampleMeta = samples.merge(meta, how = "left", on = "ID")


complex_schema_pop = tskit.MetadataSchema({
    'codec': 'json',
    'additionalProperties': False,
    'properties': {'pop': {'description': 'Lineage',
                           'type': 'string'},
                   'subpop': {'description': 'Subspecies',
                           'type': 'string'}},
    'required': ['pop'],
    'type': 'object',
})

tables = ts.dump_tables()
tables.populations.metadata_schema = complex_schema_pop

# Add data to population table
popdata_column = [{"pop": pop, "subpop":subpop} for pop,subpop in set(list(zip(sampleMeta.Pop, sampleMeta.SubPop)))]
encoded_popdata_column = [
    tables.populations.metadata_schema.validate_and_encode_row(r) for r in popdata_column
]
md, md_offset = tskit.pack_bytes(encoded_popdata_column)
tables.populations.set_columns(metadata=md, metadata_offset=md_offset)


# %%
popDict = {y: x for x,y in enumerate(set(list(zip(sampleMeta.Pop, sampleMeta.SubPop))))}
ind_nodes = list([x.id for x in ts.nodes() if x.time == 0.0])
ind_pops = [popDict[(sampleMeta.Pop.to_list()[x], sampleMeta.SubPop.to_list()[x])] for x in ind_nodes]


# %%
indNodePop = tables.nodes.population
indNodePop[ind_nodes] = ind_pops
tables.nodes.population = indNodePop
ts_meta = tables.tree_sequence()

#Remove reported siblings
removeNames = pd.read_csv("/home/jobste/Singer/RemovedSibling.txt", names = ["ID"])
keepNodes = [num for num,name in enumerate(sampleMeta.ID.to_list()) if name not in removeNames.ID.to_list()]
ts_meta = ts_meta.simplify(keepNodes)


ts_meta = tables.tree_sequence()
lineages = set([x.metadata["pop"] for x in ts_meta.populations()])
subspecies = set([x.metadata["subpop"] for x in ts_meta.populations()])

nodes_by_ind = [[x] for x in ind_nodes]

#GNN
samples_listed_by_subsp = {
    pop: [node for node,nodePop in zip(ind_nodes, ind_pops) if ts_meta.population(nodePop).metadata["subpop"] == pop] for pop in subspecies
}
#samples_listed_by_subsp = [x for x in samples_listed_by_subsp if len(x) != 0]

samples_listed_by_lineage = {
    pop: [node for node,nodePop in zip(ind_nodes, ind_pops) if ts_meta.population(nodePop).metadata["pop"] == pop] for pop in lineages
}

popsPairs = [(x,y) for x in range(ts_meta.num_populations) for y in range(x, ts_meta.num_populations)]


windowsList = np.linspace(0, ts_meta.sequence_length, num = 1000)

# %%
popsTajima = pd.DataFrame()
for pop in list(samples_listed_by_subsp.keys()):
    windowTajima = ts_meta.Tajimas_D(sample_sets=samples_listed_by_subsp[pop], windows=windowsList)
    tmp = pd.DataFrame({"Breakpoint": windowsList[:-1], "TajimaD": list(windowTajima), "Pop": pop})
    popsTajima = pd.concat([popsTajima, tmp])

popsTajima.to_csv("WindowedTajima_" + str(chromosome) + ".csv", index = None)


# %%
popsDiversity = pd.DataFrame()
for pop in samples_listed_by_subsp.keys():
    windowDiversity = ts_meta.diversity(sample_sets=samples_listed_by_subsp[pop], windows=windowsList)
    tmp = pd.DataFrame({"Breakpoint": windowsList[:-1], "Diversity": list(windowDiversity), "Pop": pop})
    popsDiversity = pd.concat([popsDiversity, tmp])

popsDiversity.to_csv("WindowedDiversity_" +str(chromosome) + ".csv", index = None)


# %%
popsFst = pd.DataFrame()
for pop1, pop2 in popsPairs:
    pop1Name = list(samples_listed_by_subsp.keys())[pop1]
    pop2Name = list(samples_listed_by_subsp.keys())[pop2]
    windowFst = ts_meta.Fst(sample_sets=[samples_listed_by_subsp[pop1Name], samples_listed_by_subsp[pop2Name]], windows=windowsList)
    tmp = pd.DataFrame({"Breakpoint": windowsList[:-1], "Fst": list(windowFst), "Pop1": pop1Name, "Pop2": pop2Name})
    popsFst = pd.concat([popsFst, tmp])



popsFst.to_csv("WindowedFst_" + str(chromosome) + ".csv", index = None)


popsDivergence = pd.DataFrame()
for pop1, pop2 in popsPairs:
    pop1Name = list(samples_listed_by_subsp.keys())[pop1]
    pop2Name = list(samples_listed_by_subsp.keys())[pop2]
    windowDivergence = ts_meta.divergence(sample_sets=[samples_listed_by_subsp[pop1Name], samples_listed_by_subsp[pop2Name]], windows=windowsList)
    tmp = pd.DataFrame({"Breakpoint": windowsList[:-1], "Divergence": list(windowDivergence), "Pop1": pop1Name, "Pop2": pop2Name})
    popsDivergence = pd.concat([popsDivergence, tmp])



popsDivergence.to_csv("WindowedDivergence_" + str(chromosome) + ".csv", index = None)




