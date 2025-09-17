# %% [markdown]
# 

# %%
import tskit 
import pandas as pd
import numpy as np
import pickle
import sys
import cloudpickle

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
                           'type': 'string'},
                   'country': {'description': 'Country',
                           'type': 'string'}},
    'required': ['pop'],
    'type': 'object',
})

tables = ts.dump_tables()
tables.populations.metadata_schema = complex_schema_pop

# Add data to population table
popdata_column = [{"pop": pop, "subpop":subpop, "country": country } for pop,subpop,country in set(list(zip(meta.Pop, meta.SubPop, meta.Country)))]
encoded_popdata_column = [
    tables.populations.metadata_schema.validate_and_encode_row(r) for r in popdata_column
]
md, md_offset = tskit.pack_bytes(encoded_popdata_column)
tables.populations.set_columns(metadata=md, metadata_offset=md_offset)

popDict = {y: x for x,y in enumerate(set(list(zip(meta.Pop, meta.SubPop, meta.Country))))}
ind_nodes = list([x.id for x in ts.individuals()])
ind_pops = [popDict[(sampleMeta.Pop.to_list()[x], sampleMeta.SubPop.to_list()[x], sampleMeta.Country.to_list()[x])] for x in ind_nodes]

# %%
indNodePop = tables.nodes.population
indNodePop[ind_nodes] = ind_pops
tables.nodes.population = indNodePop
ts_meta = tables.tree_sequence()

#Remove reported siblings
removeNames = pd.read_csv("/home/jobste/Singer/RemovedSibling.txt", names = ["ID"])
keepNodes = [num for num,name in enumerate(sampleMeta.ID.to_list()) if name not in removeNames.ID.to_list()]
ts_meta = ts_meta.simplify(keepNodes)


# %%
subspecies_list = [list(ts_meta.samples(popNum)) for popNum in range(ts_meta.num_populations)]
subspecies_listName = [ts_meta.population(popNum).metadata["subpop"] for popNum in range(ts_meta.num_populations)]

pd.DataFrame({"Subspecies": subspecies_listName}).to_csv("Chr" + str(chromosome) + "_CountriesNames.csv", index = None)

# %%
coaltime = ts_meta.coalescence_time_distribution(sample_sets=subspecies_list, weight_func="pair_coalescence_events") # (0,0) (0, 1) (0, 2) (1, 1), (1, 2)
cloudpickle.dump(coaltime, open("CoalRateRaw_Country_" + str(chromosome) + ".pkl", "wbr"))

