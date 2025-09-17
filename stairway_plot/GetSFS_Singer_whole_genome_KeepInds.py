

# %%
import tskit
import pandas as pd
import sys
import numpy as np
from math import ceil

chromosomes = range(1, 17)

popType = sys.argv[1]
popName = sys.argv[2]
foldType = sys.argv[3]
if len(sys.argv) > 4:
    dmax_cutoff = sys.argv[4]


print(popType + "_" + popName)


def multimap(s):
    if len(s.mutations) == 2:
        if not np.any([m.edge == tskit.NULL for m in s.mutations]):
            return s.id
    elif len(s.mutations) > 2:
        return s.id



# Read in meta
# %%
meta = pd.read_csv("/sietch_colab/jobste/Singer/NewRun_CleanInds/TsinferMetaData_Vcf_Country_corrected.txt")

samples = pd.read_csv("/sietch_colab/jobste/Singer/NewRun_CleanInds/SingerSamples_haploid.txt", header = None, names = ["ID"])
sampleMeta = samples.merge(meta, how = "left", on = "ID")

for num,chromosome in enumerate(chromosomes):
    # %%
    ts = tskit.load("/sietch_colab/jobste/Singer/NewRun_CleanInds/Chr" + str(chromosome) + "_Singer_99.trees")

    # Remove error mutations
    remove_sites = [multimap(s) for s in ts.sites()]
    remove_sites = [x for x in remove_sites if x is not None]
    ts = ts.delete_sites(remove_sites)


    # %%
    # Add subpop as pop
    ######################################################
    if popType == "subpop":
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

        popDict = {y: x for x,y in enumerate(set(list(zip(sampleMeta.Pop, sampleMeta.SubPop))))}
        ind_nodes = list([x.id for x in ts.nodes() if x.time == 0.0])
        ind_pops = [popDict[(sampleMeta.Pop.to_list()[x], sampleMeta.SubPop.to_list()[x])] for x in ind_nodes]
    ######################################################3

    # Add lineage as pop
    ######################################################
    if popType == "pop":
        complex_schema_pop = tskit.MetadataSchema({
            'codec': 'json',
            'additionalProperties': False,
            'properties': {'pop': {'description': 'Lineage',
                                'type': 'string'}},
            'required': ['pop'],
            'type': 'object',
        })

        tables = ts.dump_tables()
        tables.populations.metadata_schema = complex_schema_pop

        # Add data to population table
        popdata_column = [{"pop": pop} for pop in set(sampleMeta.Pop)]
        encoded_popdata_column = [
            tables.populations.metadata_schema.validate_and_encode_row(r) for r in popdata_column
        ]
        md, md_offset = tskit.pack_bytes(encoded_popdata_column)
        tables.populations.set_columns(metadata=md, metadata_offset=md_offset)

        popDict = {y: x for x,y in enumerate(set(sampleMeta.Pop))}
        ind_nodes = list([x.id for x in ts.nodes() if x.time == 0.0])
        ind_pops = [popDict[sampleMeta.Pop.to_list()[x]] for x in ind_nodes]
    ######################################################3


    indNodePop = tables.nodes.population
    indNodePop[ind_nodes] = ind_pops
    tables.nodes.population = indNodePop
    ts_meta = tables.tree_sequence()

    if len(sys.argv) > 4:
        print("Keeping only clean individuals.")
        # Simplify to the clean samples
        keep = pd.read_csv("Keep_" + popName + "_dmax_" + str(dmax_cutoff) + "_highFreq.csv", header=None)[0].to_list()
        keep_nodes = sampleMeta[sampleMeta.ID.isin(keep)].index
        print(len(keep_nodes))
    else: 
        print("Keeping all individuals.")
        popNum = [x.id for x in ts_meta.populations() if x.metadata[popType] == popName][0]
        keep_nodes = list(ts_meta.samples(popNum))
        print(len(keep_nodes))

    if foldType == "unfold":
        if num == 0:
            sfs = ts_meta.allele_frequency_spectrum(sample_sets = [keep_nodes], polarised = True, span_normalise = False)
        else:
            tmpSfs = ts_meta.allele_frequency_spectrum(sample_sets = [keep_nodes], polarised = True, span_normalise = False)
            sfs = sfs + tmpSfs


    elif foldType == "fold":
        if num == 0:
            sfs = ts_meta.allele_frequency_spectrum(sample_sets = [keep_nodes], polarised = False, span_normalise = False)
        else:
            tmpSfs = ts_meta.allele_frequency_spectrum(sample_sets = [keep_nodes], polarised = False, span_normalise = False)
            sfs = sfs + tmpSfs

if foldType == "unfold":
    pd.DataFrame({"Freq": sfs[1:-1]}).to_csv("Genome_pop" + str(popName) + "_unfold_count.sfs", index = None, header=None)
elif foldType == "fold":
    pd.DataFrame({"Freq": sfs[1:ceil(len(keep_nodes)/2)]}).to_csv("Genome_pop" + str(popName) + "_unfold_count.sfs", index = None, header=None)



# %
# %%
noSamples = len(keep_nodes)
pd.DataFrame({"NoSamples": [noSamples]}).to_csv("NoSamples_Pop" + popName + ".csv", index = None, header = None)
pd.DataFrame({"NoSamples": [round((noSamples - 2) / 4)]}).to_csv("NoSamples_Pop1" + popName + ".csv", index = None, header = None)
pd.DataFrame({"NoSamples": [round((noSamples - 2) / 2)]}).to_csv("NoSamples_Pop2" + popName + ".csv", index = None, header = None)
pd.DataFrame({"NoSamples": [round((noSamples - 2) * 0.75)]}).to_csv("NoSamples_Pop3" + popName + ".csv", index = None, header = None)
pd.DataFrame({"NoSamples": [round((noSamples - 2) )]}).to_csv("NoSamples_Pop4" + popName + ".csv", index = None, header = None)
