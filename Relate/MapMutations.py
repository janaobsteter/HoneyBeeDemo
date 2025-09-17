# %%
import tskit
import pandas as pd
import allel
import numpy as np
from collections import Counter
from collections import defaultdict
import matplotlib.pyplot as plt
import sys
from math import floor

chromosome = sys.argv[1]
chunkTrees = int(sys.argv[2])

# %%
# Read the tree sequence
chrTree = tskit.load("/home/jobste/Relate/All_snps_with_outgroup_info_newAA/Chr" + str(chromosome) + ".trees")
tables = chrTree.dump_tables()

# %%
# Read vcf
vcf = allel.read_vcf("../Data/Chr" + str(chromosome) + "_imputed_AllInfo.vcf.gz")
# Get positions and genotypes
pos = vcf['variants/POS']
gt = vcf['calldata/GT']
ref = vcf['variants/REF']
alt = vcf['variants/ALT'][:, 0]




# Loop over trees
noHits = defaultdict()

chunkSize = floor(chrTree.num_trees / 100)
chunkStart = chunkTrees * chunkSize
chunkEnd = (chunkTrees + 1) * chunkSize
if int(chunkTrees) != 99:
    selectedTrees = range(chunkStart, chunkEnd)
else:
    selectedTrees = range(chunkStart, chrTree.num_trees+1)

for treeIndex in selectedTrees:
    oneTree = chrTree.at_index(treeIndex)

    intervalBoolean = (vcfPos > oneTree.interval.left) & (vcfPos < oneTree.interval.right)
    intervalPositions = vcfPos[intervalBoolean]
    missingBoolean = list(map(lambda x: x not in treePos, intervalPositions))
    missingGenotypes = gt[intervalBoolean][:, :, 0][missingBoolean]
    missingPositions = pos[intervalBoolean][missingBoolean]
    missingAlleles = list(zip(ref[intervalBoolean][missingBoolean], alt[intervalBoolean][missingBoolean]))
    for seq, site in enumerate(missingPositions):
        alleles = list(set(missingAlleles[seq]))
        anc1, parsMut1 = oneTree.map_mutations(genotypes=missingGenotypes[seq], alleles = missingAlleles[seq], ancestral_state = alleles[0])
        anc2, parsMut2 = oneTree.map_mutations(genotypes=missingGenotypes[seq], alleles = missingAlleles[seq], ancestral_state = alleles[1])
        
        if len(parsMut1) > len(parsMut2):
            ancestral = anc2
            parsMut = parsMut2
            polarisation = "Alt"
        elif len(parsMut2) > len(parsMut1):
            ancestral = anc1
            parsMut = parsMut1
            polarisation = "Ref"
        else:
            continue
    
        site_id = tables.sites.add_row(site, ancestral_state=ancestral)
        parent_offset = len(tables.mutations)
        for mut in parsMut:
            parent = mut.parent
            if parent != tskit.NULL:
                parent += parent_offset
            mut_id = tables.mutations.add_row(
                site_id, node=mut.node, parent=parent, derived_state=mut.derived_state)


        noHits[site] = {"TreeIndex": treeIndex, 
                        "NoHits": len(parsMut),
                        "Polarisation": polarisation,
                        "AF": sum(missingGenotypes[seq]) / len(missingGenotypes[seq]),
                        "NumSamples": {parsimonious_mut.node: oneTree.num_samples(parsimonious_mut.node) for parsimonious_mut in parsMut}}


# Create new tree from the added mutations
tables.sort()        
new_ts = tables.tree_sequence()
new_ts.dump("/home/jobste/Relate/All_snps_with_outgroup_info_newAA/Chr" + str(chromosome) + "_chunk" + str(chunkTrees) + "_remapped.trees")
# Dataframe of number of hits per site
pd.DataFrame.from_records([(x, y['NoHits'], y['Polarisation'], y["AF"], y['NumSamples']) for x, y in noHits.items()], 
                          columns = ['Site', 'NoHits', 'Polarisation', 'AF', 'NumSamples']).\
                            to_csv("/home/jobste/Relate/All_snps_with_outgroup_info_newAA/MappedMutations" + str(chromosome) + "_chunk" + str(chunkTrees) + ".csv", index = False)




