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

tree = tskit.load("All_snps_with_outgroup_info_newAA/Chr" + str(chromosome) + ".trees")
treePos = [x.position for x in tree.sites()]
chrSample = pd.read_table("All_snps_with_outgroup_info_newAA/Chr" + str(chromosome) + "_relate.sample")


# %%
meta = pd.read_csv("TsinferMetaData_Vcf.txt")
meta.head()

sampleMeta = chrSample.merge(meta, how = "left", left_on = "ID_1", right_on = "ID")
sampleMeta = sampleMeta.iloc[1: , :]

complex_schema = tskit.MetadataSchema({
    'codec': 'json',
    'additionalProperties': False,
    'properties': {'name': {'description': 'ENA accession number',
                                 'type': 'string'},
                   'pop': {'description': 'Lineage',
                           'type': 'string'},
                    'subpop': {'description': 'Subspecies',
                           'type': 'string'}},
    'required': ['name', 'pop', 'subpop'],
    'type': 'object',
})


# Add meta data
tables = tree.dump_tables()
tables.individuals.metadata_schema = complex_schema



# Add metadata
metadata_column = [{"name": list(sampleMeta.ID_1)[ind], "pop": list(sampleMeta.Pop)[ind], "subpop": list(sampleMeta.SubPop)[ind]} for ind in range(tree.num_individuals)]
encoded_metadata_column = [
    tables.individuals.metadata_schema.validate_and_encode_row(r) for r in metadata_column
]
md, md_offset = tskit.pack_bytes(encoded_metadata_column)
tables.individuals.set_columns(flags=range(tree.num_individuals), metadata=md, metadata_offset=md_offset)

tree_meta = tables.tree_sequence()

# %%
sample_nodes = [tree_meta.node(n) for n in tree_meta.samples()]
sample_ids = [n.id for n in sample_nodes]
sample_names = [
    tree_meta.individual(n).metadata["name"]
    for n in sample_ids
]
sample_pops = [
    tree_meta.individual(n).metadata["pop"]
    for n in sample_ids
]
sample_subpops = [
    tree_meta.individual(n).metadata["subpop"]
    for n in sample_ids
]

pops = dict()
for (node, pop) in zip(sample_ids, sample_pops):
    if not pop in pops.keys():
        pops[pop] = [node]
    else:
        pops[pop].append(node)
numPops = len(pops.keys())

subpops = dict()
for (node, pop) in zip(sample_ids, sample_subpops):
    if not pop in subpops.keys():
        subpops[pop] = [node]
    else:
        subpops[pop].append(node)
numSubpops = len(subpops.keys())

popsPairs = [(x, y) for x in range(numPops) for y in range(numPops) if x < y]
subpopsPairs = [(x, y) for x in range(numSubpops) for y in range(numSubpops) if x < y]


# %%
tree.samples()
windowsList = np.linspace(0, tree.sequence_length, num = 1000)

# %%
popsTajima = pd.DataFrame()
for pop in pops.keys():
    windowTajima = tree.Tajimas_D(sample_sets=pops[pop], windows=windowsList)
    tmp = pd.DataFrame({"Breakpoint": windowsList[:-1], "TajimaD": list(windowTajima), "Pop": pop})
    popsTajima = pd.concat([popsTajima, tmp])

popsTajima.to_csv("WindowedTajima_" + str(chromosome) + ".csv", index = None)


# %%
popsDiversity = pd.DataFrame()
for pop in pops.keys():
    windowDiversity = tree.diversity(sample_sets=pops[pop], windows=windowsList)
    tmp = pd.DataFrame({"Breakpoint": windowsList[:-1], "Diversity": list(windowDiversity), "Pop": pop})
    popsDiversity = pd.concat([popsDiversity, tmp])

popsDiversity.to_csv("WindowedDiversity_" +str(chromosome) + ".csv", index = None)


# %%
popsFst = pd.DataFrame()
for pop1, pop2 in popsPairs:
    pop1Name = list(pops.keys())[pop1]
    pop2Name = list(pops.keys())[pop2]
    windowFst = tree.Fst(sample_sets=[pops[pop1Name], pops[pop2Name]], windows=windowsList)
    tmp = pd.DataFrame({"Breakpoint": windowsList[:-1], "Fst": list(windowFst), "Pop1": pop1Name, "Pop2": pop2Name})
    popsFst = pd.concat([popsFst, tmp])



popsFst.to_csv("WindowedFst_" + str(chromosome) + ".csv", index = None)


popsDivergence = pd.DataFrame()
for pop1, pop2 in popsPairs:
    pop1Name = list(pops.keys())[pop1]
    pop2Name = list(pops.keys())[pop2]
    windowDivergence = tree.divergence(sample_sets=[pops[pop1Name], pops[pop2Name]], windows=windowsList)
    tmp = pd.DataFrame({"Breakpoint": windowsList[:-1], "Divergence": list(windowFst), "Pop1": pop1Name, "Pop2": pop2Name})
    popsDivergence = pd.concat([popsDivergence, tmp])



popsDivergence.to_csv("WindowedDivergence_" + str(chromosome) + ".csv", index = None)



