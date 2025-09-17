# %% [markdown]
# 

# %%
import tskit 
import pandas as pd
import numpy as np
import pickle
import sys

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




# %%
coaltime_NoSamples = ts_meta.coalescence_time_distribution(weight_func="pair_coalescence_events")

# %%
time_windows = coaltime_NoSamples.quantile(np.linspace(0, 1, 25))

output = open('CoalRates_Subspecies_' + str(chromosome) + '_TimeWindow.pkl', 'wb')
pickle.dump(time_windows, output)
output.close()