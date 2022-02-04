#!/usr/bin/env python
# coding: utf-8

import cyvcf2
import tsinfer
import pandas as pd
import json
import numpy as np
import sys
import os

chromosome = sys.argv[1]
vcfFile = sys.argv[2]
metaFile = sys.argv[3]
ancestralFile = sys.argv[4]
vcfFileName = vcfFile.split("/")[-1].strip(".vcf")
sampleFile = os.getcwd() + "/" + str(vcfFileName) + ".samples"

print("Chromosome is " + str(chromosome))

#######################################################################
# Define the functions to read in the vcf
#######################################################################
# -- Haploid --#
def add_haploid_individuals(vcf, samples, populations):
    """
    The function to store the information about which nodes pertain to which sample
    Could be modified in case of a different ploidy
    :param vcf:
    :param samples: the sample file to hold tskit.SampleData
    :param populations: population of the samples
    :return:
    """
    for name, population in zip(vcf.samples, populations):
        samples.add_individual(ploidy=1, metadata={"name": name}, population=population)




def add_haploid_sites(vcf, samples, ancDF):
    """
    Read the sites in the vcf and add them to the samples object, reordering the
    alleles to put the ancestral allele first, if it is available. Read the ancestral allele from the ancestral DF.
    """
    pos = 0

    for variant in vcf:  # Loop over variants, each assumed at a unique site
        if pos == variant.POS:
            raise ValueError("Duplicate positions for variant at position", pos)
        else:
            pos = variant.POS
        #if any([not phased for _, _, phased in variant.genotypes]):
        #    raise ValueError("Unphased genotypes for variant at position", pos)
        alleles = [variant.REF] + variant.ALT
        try:
            ancestral = list(ancDF.AncAl[ancDF.SnpPos == str(chromosome) + "_" + str (pos)])[0]
        except:
            continue #ancestral = variant.INFO.get("AA", variant.REF)        # Ancestral state must be first in the allele list.
        ordered_alleles = [ancestral] + list(set(alleles) - {ancestral})
        allele_index = {
            old_index: ordered_alleles.index(allele)
            for old_index, allele in enumerate(alleles)
        }
        allele_index[-1] = -1
        # Map original allele indexes to their indexes in the new alleles list.
        genotypes = [
            allele_index[old_index]
            for row in variant.genotypes
            for old_index in row[0:1]
        ]
        samples.add_site(pos, genotypes=genotypes, alleles=ordered_alleles)

# -- Diploid --#
def add_diploid_individuals(vcf, samples, populations):
    """
    The function to store the information about which nodes pertain to which sample
    Could be modified in case of a different ploidy
    :param vcf: Vcf file of your samples
    :param samples:
    :param populations:
    :return:
    """
    for name, population in zip(vcf.samples, populations):
        samples.add_individual(ploidy=2, metadata={"name": name}, population=population)

def add_diploid_sites(vcf, samples, ancDF):
    """
    Read the sites in the vcf and add them to the samples object, reordering the
    alleles to put the ancestral allele first, if it is available.
    """
    pos = 0
    for variant in vcf:  # Loop over variants, each assumed at a unique site
        if pos == variant.POS:
            raise ValueError("Duplicate positions for variant at position", pos)
        else:
            pos = variant.POS
        if any([not phased for _, _, phased in variant.genotypes]):
            raise ValueError("Unphased genotypes for variant at position", pos)
        alleles = [variant.REF] + variant.ALT
        try:
            ancestral = list(ancDF.AncAl[ancDF.SnpPos == "carnica.LG" + str(chromosome) + "_" + str (pos)])[0]
        except:
            ancestral = variant.INFO.get("AA", variant.REF)
        # Ancestral state must be first in the allele list.
        ordered_alleles = [ancestral] + list(set(alleles) - {ancestral})
        allele_index = {
            old_index: ordered_alleles.index(allele)
            for old_index, allele in enumerate(alleles)
        }
        allele_index[-1] = -1
        # Map original allele indexes to their indexes in the new alleles list.
        genotypes = [
            allele_index[old_index]
            for row in variant.genotypes
            for old_index in row[0:2]
        ]
        samples.add_site(pos, genotypes=genotypes, alleles=ordered_alleles)


# This function checks whether there is only one contig in the vcf - this may not be the case!
def chromosome_length(vcf, chromosome):
    """
    This is a modified function to read in the sequence length
    :param vcf: vcffile
    :return: the length of the sequence (genome or chromosme - depending what the vcf is for)
    """
    #assert len(vcf.seqlens) == 1
    return vcf.seqlens[int(chromosome)-1]


def add_populations(vcf, samples):
    """
    Add tsinfer Population objects and returns a list of IDs corresponding to the VCF samples.
    """
    # In this VCF, the first letter of the sample name refers to the population
    samples_first_letter = [sample_name[0] for sample_name in vcf.samples]
    pop_lookup = {}
    pop_lookup["8"] = samples.add_population(metadata={"country": "Norway"})
    pop_lookup["F"] = samples.add_population(metadata={"country": "France"})
    print(pop_lookup)
    print(samples.metadata)
    return [pop_lookup[first_letter] for first_letter in samples_first_letter]

def add_populations_drones(vcf, samples, metaData):
    """
    This is a modified add populations function for drones from David Wragg data
    Add tsinfer Population objects for drone subspecies, stored in metaPop object, and return a list od IDs correposning to the VCF samples
    Input: vcf = cyvcf2 VCF object, samples is tsinfer SampleData and metaData is a pandas dataframe with id in ID column and populations in Type column
    Return: A list of population indexes
    """
    pop_lookup = {}
    pop_lookup["carnica"] = samples.add_population(metadata={"subspecies": "carnica", "lineage": "C"})
    pop_lookup["mellifera"] = samples.add_population(metadata={"subspecies": "mellifera", "lineage": "M"})
    pop_lookup["scutellata"] = samples.add_population(metadata={"subspecies": "scutellata", "lineage": "A"})
    pop_lookup["caucasica"] = samples.add_population(metadata={"subspecies": "caucasica", "lineage": "O"})
    pop_lookup["unicolor"] = samples.add_population(metadata={"subspecies": "unicolor", "lineage": "A"})
    pop_lookup["hybrid"] = samples.add_population(metadata={"subspecies": "hybrid", "lineage": "H"})
    pop_lookup["capensis"] = samples.add_population(metadata={"subspecies": "capensis", "lineage": "A"})
    pop_lookup["ligustica"] = samples.add_population(metadata={"subspecies": "ligustica", "lineage": "C"})
    return [pop_lookup[list(metaData.Type[metaData.ID == x])[0]] for x in vcf.samples]


#######################################################################
#######################################################################

#######################################################################
# Read in the information and prepare the data
#######################################################################
print("Read in the vcf file")
# The files with the meta data about the samples, including the subspecie (Type)
meta = pd.read_csv(metaFile)
meta = meta[['ID', 'Type']]
ancAl = pd.read_csv(ancestralFile, sep=",", header=None)
ancAl.columns = ["SnpPos", "AncAl"]
#meta.tail(50)
#print(type(meta.Type[meta.ID == "scu_d_2_3"]))

# Create a population (subspecie) list for the samples in the VCF
vcfD = cyvcf2.VCF(vcfFile, strict_gt=True)
metaPop = [list(meta.Type[meta.ID == x])[0] if x in list(meta.ID) else "hybrid" for x in vcfD.samples]
metaPopDF = pd.DataFrame({'ID': [x for x in vcfD.samples], 'Type': metaPop})

# Create tsinfer.SampleData object
# Modify the length of the sequence - it can either be read from the vcf files (if there is only one contig), otherwise adjust manually
#-- For haploid data --#
vcfD = cyvcf2.VCF(vcfFile, strict_gt=True)

#with tsinfer.SampleData(
#        path=sampleFile, 
#	sequence_length=chromosome_length(vcfD, chromosome),
#	num_flush_threads=2,
#	max_file_size=2**30,
#) as samples:
#    populations = add_populations_drones(vcfD, samples, metaPopDF)
#    print("populations determined")
#    add_haploid_individuals(vcfD, samples, populations)
#    print("Inds added")
#    add_haploid_sites(vcfD, samples, ancAl)

#print(
#    "Sample file created for {} samples ".format(samples.num_samples)
#    + "({} individuals) ".format(samples.num_individuals)
#    + "with {} variable sites.".format(samples.num_sites)
#)

#-- For diploid data --#
#with tsinfer.SampleData(
#        path="DronesD.samples", sequence_length=7238532,
#) as samples:
#    populations = add_populations_drones(vcfD, samples, meta)
#    add_diploid_individuals(vcfD, samples, populations)
#    add_diploid_sites(vcfD, samples)

#print(
#    "Sample file created for {} samples ".format(samples.num_samples)
#    + "({} individuals) ".format(samples.num_individuals)
#    + "with {} variable sites.".format(samples.num_sites)
#)


#######################################################################
# Do the inference and write the outputs
#######################################################################

# Do the inference on the 10 SNPs
samples = tsinfer.load("Chr" + str(chromosome) + ".samples")

drone_ts = tsinfer.infer(samples)
print(
    "Inferred tree sequence `{}`: {} trees over {} Mb".format(
        "drone_ts", drone_ts.num_trees, drone_ts.sequence_length / 1e6
    )
)

# Check the metadata
for sample_node_id in drone_ts.samples():
    individual_id = drone_ts.node(sample_node_id).individual
    population_id = drone_ts.node(sample_node_id).population
    print(
        "Node",
        sample_node_id,
        "labels genome sampled from",
        json.loads(drone_ts.individual(individual_id).metadata),
        "in",
        json.loads(drone_ts.population(population_id).metadata)["subspecies"],
    )

drone_ts.dump("Chr" + str(chromosome) + ".trees")


print("Draw a tree")
# Set the colours of the subspecies for plotting
colours = {"hybrid": "grey", "carnica": "yellow", "ligustica": "orange", "mellifera": "red", "caucasica": "blue", "scutellata": "black", "unicolor": "blue", "capensis": "purple", }
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

# Draw a tree at 1e6, coloured by subspecies
#tree = drone_ts.at(1e6)
#tree.draw(
#    path="tree_at_1Mb" + str(chromosome) + ".svg",
#    height=700,
#    width=1200,
#    node_labels=individual_for_node,
#    node_colours=colours_for_node
#)

# Draw a tree at 1e6, coloured by lineages
#tree.draw(
#    path="tree_at_1Mb_lineage" + str(chromosome) + ".svg",
#    height=700,
#    width=1200,
#    node_labels=individual_for_node,
#    node_colours=colours_for_nodeL
#)


# In[198]:


#Obtain node info - but just for samples!
# There is also intermediate nodes!!! (total #nodes > #nodes for samples)
sample_nodes = [drone_ts.node(n) for n in drone_ts.samples()]

#Get samples ids
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



# In[141]:


print("Do the statistics.")
# Obtain Fst and write a table
Fst = drone_ts.Fst([pops["A"], pops["C"], pops["M"], pops["O"]], indexes=[(0,1), (0,2), (0, 3), (1, 2), (1, 3), (2,3)])
#print(Fst)
fstDF = pd.DataFrame(Fst, columns = ["Fst"])
fstDF.loc[:, "Group1"] = ["A", "A", "A", "C", "C", "M"]
fstDF.loc[:, "Group2"] = ["C", "M", "O", "M", "O", "O"]
fstDF.to_csv("FST" + str(chromosome) + ".csv", index=None)


# In[144]:


# Obtain Tajima's D and write a table
tajimaD = drone_ts.Tajimas_D([pops["A"], pops["C"], pops["M"], pops["O"], pops["H"]])
#print(tajimaD)
tajimaDF = pd.DataFrame(tajimaD, columns = ["TajimasD"])
tajimaDF.loc[:, "Group"] = ["A", "C", "M", "O", "H"]
tajimaDF.to_csv("TajimasD" + str(chromosome) + ".csv", index=None)


# In[145]:


# Divergence
divergence = drone_ts.divergence([pops["A"], pops["C"], pops["M"], pops["O"]], indexes=[(0,1), (0,2), (0, 3), (1, 2), (1, 3), (2,3)])
#print(divergence)
divergenceDF = pd.DataFrame(divergence, columns = ["Divergence"])
divergenceDF.loc[:, "Group1"] = ["A", "A", "A", "C", "C", "M"]
divergenceDF.loc[:, "Group2"] = ["C", "M", "O", "M", "O", "O"]
divergenceDF.to_csv("Divergence" + str(chromosome) + ".csv", index=None)



# Diversity and write a table
diversity = drone_ts.diversity([pops["A"], pops["C"], pops["M"], pops["O"], pops["H"]])
diversityDF = pd.DataFrame(diversity, columns = ["Diversity"])
diversityDF.loc[:, "Group"] = ["A", "C", "M", "O", "H"]
diversityDF.to_csv("Diversity" + str(chromosome) + ".csv", index=None)


# Allele frequency spectrum
np.savetxt("AFS_Alineage" + str(chromosome) + ".csv", drone_ts.allele_frequency_spectrum([pops["A"]]))
np.savetxt("AFS_Mlineage" + str(chromosome) + ".csv", drone_ts.allele_frequency_spectrum([pops["M"]]))
np.savetxt("AFS_Clineage" + str(chromosome) + ".csv", drone_ts.allele_frequency_spectrum([pops["C"]]))
np.savetxt("AFS_Olineage" + str(chromosome) + ".csv", drone_ts.allele_frequency_spectrum([pops["O"]]))


# In[228]:


# f2, f3, f4 statistics
f2 = drone_ts.f2([pops["A"], pops["C"], pops["M"], pops["O"]], indexes=[(0,1), (0,2), (0, 3), (1, 2), (1, 3), (2,3)])
#print(f2)
f2DF = pd.DataFrame(f2, columns = ["f2"])
f2DF.loc[:, "Group1"] = ["A", "A", "A", "C", "C", "M"]
f2DF.loc[:, "Group2"] = ["C", "M", "O", "M", "O", "O"]
f2DF.to_csv("f2" + str(chromosome) + ".csv", index=None)

# f4
f4 = drone_ts.f4([pops["A"], pops["C"], pops["M"], pops["O"]])
print("This is f4 statistics: ", str(f4))



# GNN
gnn = pd.DataFrame(drone_ts.genealogical_nearest_neighbours(pops["H"], sample_sets=[pops["A"], pops["C"], pops["M"], pops["O"]]), columns = ["A", "C","M", "O"])

# Add sample ids and names to the data frame
#gnn.loc[:, "Id"] = pops["H"] #The two nodes per individual will be the same since they are homozygous diploids
namesH = [name for (name, ID) in zip(sample_names, sample_ids) if ID in pops['H']]
gnn.loc[:, "Names"] = namesH
gnn = gnn.drop_duplicates()

gnn.to_csv("GNN_drones" + str(chromosome) + ".csv", index=None)

# Genetic relatedness
grel = drone_ts.genetic_relatedness([pops["A"], pops["C"], pops["M"], pops["O"]], indexes=[(0,1), (0,2), (0, 3), (1, 2), (1, 3), (2,3)])
#print(grel)
grelDF = pd.DataFrame(grel, columns = ["GeneticRel"])
grelDF.loc[:, "Group1"] = ["A", "A", "A", "C", "C", "M"]
grelDF.loc[:, "Group2"] = ["C", "M", "O", "M", "O", "O"]
grelDF.to_csv("GeneticRelatedness" + str(chromosome) + ".csv", index=None)


# In[221]:


# Mean descendants of ALL nodes!!! Also intermediate
meanDesc = pd.DataFrame(drone_ts.mean_descendants([pops["A"], pops["C"], pops["M"], pops["O"]]), columns = ["A", "C", "M", "O"])
meanDesc.loc[:, "NodeID"] = [n.id for n in drone_ts.nodes()]
meanDesc.loc[:, "NodeName"] = None

# Add names for sample nodes
for (ID, name) in zip(sample_nodes, sample_names):
    meanDesc.at[meanDesc.NodeID == ID.id, "NodeName"] = name


meanDesc.to_csv("MeanDescendants" + str(chromosome) + ".csv", index=None)



# Write the tables with mutations, nodes, sites, populations and individuals
drone_ts.dump_text(mutations=open("Mutations" + str(chromosome) + ".txt", "w"), nodes = open("Nodes" + str(chromosome) + ".txt", "w"), sites = open("Sites" + str(chromosome) + ".txt", "w"), edges = open("Edges" + str(chromosome) + ".txt", "w"), individuals = open("Individuals" + str(chromosome) + ".txt", "w"), populations = open("Populations" + str(chromosome) + ".txt", "w", encoding="utf-8"))

