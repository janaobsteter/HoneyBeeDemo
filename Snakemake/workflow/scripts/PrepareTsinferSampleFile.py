#!/usr/bin/env python
# coding: utf-8

import cyvcf2
import tsinfer
import pandas as pd
import json
import numpy as np
import sys
import os

chromosome = snakemake.wildcards['chromosome']
vcfFile = snakemake.input[0]
meta = snakemake.input[1]
# TODO: Prepare a file with major alleles
major = snakemake.input[2]
pops = snakemake.config['populations']
sampleFile = snakemake.output[0]


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

def add_haploid_sites(vcf, samples):
    """
    Read the sites in the vcf and add them to the samples object, reordering the
    alleles to put the ancestral allele first, if it is available. Read the ancestral
    allele from the the vcf AA field
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
        # try:
        #     ancestral = list(ancDF.AncAl[ancDF.SnpPos == str(chromosome) + "_" + str (pos)])[0]
        # except:
        #     continue
        ancestral = variant.INFO.get("AA", variant.REF)
        if ancestral == ".":
            varPos = str(chromosome) + "_" + str(var.POS)
            ancestral = list(major.Anc[major.Pos == varPos])[0]

        try:
            ancestral_allele = alleles.index(ancestral)
            #print(pos, ancestral, ancestral_allele)
        except:
            ancestral_allele = MISSING_DATA

        genotypes = [g for row in variant.genotypes for g in row[0:2]]
        #print('POS: ', pos, ' ANC: ', ancestral, ' ALLE: ', alleles)
        sample_data.add_site(position=pos,
                         genotypes=genotypes,
                         alleles=alleles,
                         ancestral_allele=ancestral_allele

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

# def add_diploid_sites(vcf, samples, ancDF):
#     """
#     Read the sites in the vcf and add them to the samples object, reordering the
#     alleles to put the ancestral allele first, if it is available.
#     """
#     pos = 0
#     for variant in vcf:  # Loop over variants, each assumed at a unique site
#         if pos == variant.POS:
#             raise ValueError("Duplicate positions for variant at position", pos)
#         else:
#             pos = variant.POS
#         if any([not phased for _, _, phased in variant.genotypes]):
#             raise ValueError("Unphased genotypes for variant at position", pos)
#         alleles = [variant.REF] + variant.ALT
#         try:
#             ancestral = list(ancDF.AncAl[ancDF.SnpPos == "carnica.LG" + str(chromosome) + "_" + str (pos)])[0]
#         except:
#             ancestral = variant.INFO.get("AA", variant.REF)
#         # Ancestral state must be first in the allele list.
#         ordered_alleles = [ancestral] + list(set(alleles) - {ancestral})
#         allele_index = {
#             old_index: ordered_alleles.index(allele)
#             for old_index, allele in enumerate(alleles)
#         }
#         allele_index[-1] = -1
#         # Map original allele indexes to their indexes in the new alleles list.
#         genotypes = [
#             allele_index[old_index]
#             for row in variant.genotypes
#             for old_index in row[0:2]
#         ]
#         samples.add_site(pos, genotypes=genotypes, alleles=ordered_alleles)


# This function checks whether there is only one contig in the vcf - this may not be the case!
def chromosome_length(vcf, chromosome):
    """
    This is a modified function to read in the sequence length
    :param vcf: vcffile
    :return: the length of the sequence (genome or chromosme - depending what the vcf is for)
    """
    #assert len(vcf.seqlens) == 1
    return vcf.seqlens[int(chromosome)-1]


# def add_populations(vcf, samples):
#     """
#     Add tsinfer Population objects and returns a list of IDs corresponding to the VCF samples.
#     """
#     # In this VCF, the first letter of the sample name refers to the population
#     samples_first_letter = [sample_name[0] for sample_name in vcf.samples]
#     pop_lookup = {}
#     pop_lookup["8"] = samples.add_population(metadata={"country": "Norway"})
#     pop_lookup["F"] = samples.add_population(metadata={"country": "France"})
#     print(pop_lookup)
#     print(samples.metadata)
#     return [pop_lookup[first_letter] for first_letter in samples_first_letter]

def add_populations_drones(vcf, samples, pops):
    """
    This is a modified add populations function for drones from David Wragg data
    Add tsinfer Population objects for drone subspecies, stored in metaPop object, and return a list od IDs correposning to the VCF samples
    Input: vcf = cyvcf2 VCF object, samples is tsinfer SampleData and metaData is a pandas dataframe with id in ID column and populations in Type column
    pops = dictinary holding name and additional information about populations
    Return: A list of population indexes
    """
    pop_lookup = {}
    for pop, lineage in pops.items():
        pop_lookup[pop] = samples.add_population(metadata={"subspecies": pop, "lineage": lineage})
    return [pop_lookup[list(metaData.Type[metaData.ID == x])[0]] for x in vcf.samples]

#######################################################################
# Read in the information and prepare the data
#######################################################################
print("Read in the vcf file")
# The files with the meta data about the samples, including the subspecie (Type)
meta = pd.read_csv(metaFile)
meta = meta[['ID', 'Type']]

# Create a population (subspecie) list for the samples in the VCF
vcfD = cyvcf2.VCF(vcfFile, strict_gt=True)
metaPop = [list(meta.Type[meta.ID == x])[0] if x in list(meta.ID) else "hybrid" for x in vcfD.samples]
metaPopDF = pd.DataFrame({'ID': [x for x in vcfD.samples], 'Type': metaPop})

# Create tsinfer.SampleData object
# Modify the length of the sequence - it can either be read from the vcf files (if there is only one contig), otherwise adjust manually
#-- For haploid data --#
vcfD = cyvcf2.VCF(vcfFile, strict_gt=True)

# Create samples for haploid data
with tsinfer.SampleData(
    path=sampleFile,
	sequence_length=chromosome_length(vcfD, chromosome),
	num_flush_threads=2,
	max_file_size=2**30,
) as samples:
   populations = add_populations_drones(vcfD, samples, metaPopDF)
   print("populations determined")
   add_haploid_individuals(vcfD, samples, populations)
   print("Inds added")
   add_haploid_sites(vcfD, samples, ancAl)

print(
   "Sample file created for {} samples ".format(samples.num_samples)
   + "({} individuals) ".format(samples.num_individuals)
   + "with {} variable sites.".format(samples.num_sites)
)

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
