#!/usr/bin/env python
# coding: utf-8

import tsdate
import tsinfer
from tskit import MISSING_DATA

import sys
import cyvcf2
import pandas as pd
import tqdm as tqdm
from math import nan

##############################################################################################
# Functions and other objects necessary to create samples ----------------------------------
##############################################################################################
def add_diploid_sites(vcf, sample_data, ancestral_state_file):
    """
    Read the sites in the vcf and add them to the samples object, reordering the
    alleles to put the ancestral allele first, if it is available.
    
    """
    pos = 0 
    progressbar = tqdm.tqdm(total=sample_data.sequence_length, desc="Read VCF", unit='bp')

    for variant in vcf:  # Loop over variants, each assumed at a unique site
        progressbar.update(variant.POS - pos)
        
        # some checks to make sure it can add sample to tsinfer
        if pos == variant.POS:
            raise ValueError("Duplicate positions for variant at position", pos)
        else:
            pos = variant.POS
        if any([not phased for _, _, phased in variant.genotypes]): 
            raise ValueError("Unphased genotypes for variant at position", pos)
        
        # Calculate allele frequency
        try:
            freq = variant.INFO.get('AC')/variant.INFO.get('AN')
        except:
            freq = nan 
        
        ref = [list(vr) for vr in variant.REF][0]
        alt = [list(va) for va in variant.ALT][0]
        alleles = ref + alt

        # Conditions to ignore sites (tsinfer don't deal with multi-alilic sites)
        if len(alleles) > 2:
            continue
        if alleles[0] == alleles[1]:
            continue
        
        # Get ancestral state from file. If position not listed use major allele.
        try:
            ancestral = list(ancestral_state_file.loc[ancestral_state_file.pos.isin([pos]), 'ances'])[0]

        except:
            if freq > 0.8:
                ancestral = ref

            elif freq < 0.2:
                ancestral = alt
            else:
                ancestral = nan
        
        # Get ancestral allele index
        try:
            ancestral_allele = alleles.index(ancestral[0])
            #print(pos, ancestral, ancestral_allele)
        except:
            ancestral_allele = MISSING_DATA
    
        genotypes = [g for row in variant.genotypes for g in row[0:2]]
        #print('POS: ', pos, ' ANC: ', ancestral, ' ALLE: ', alleles)
        sample_data.add_site(position=pos, 
                         genotypes=genotypes, 
                         alleles=alleles, 
                         ancestral_allele=ancestral_allele
                        )

def add_populations_bulls(vcf, samples, metadata):
    pop_lookup = {}
    metaPop = metadata[['Breed', 'Species']].drop_duplicates().dropna(axis = 0, how = 'all')
    sample_breed = [list(metadata.Breed[metadata.ID == x]) for x in vcf.samples]

    for breed, species in zip(metaPop.Breed, metaPop.Species):
        pop_lookup[breed] = samples.add_population(metadata={'breed':breed, 'species':species})

    return [pop_lookup[breed[0]] for breed in sample_breed]

def add_diploid_individuals(vcf, samples, populations):
    '''''
    Add metadata to sample. Because there's no pedigree, all samples are considered 'modern' except from Auroch
    which correspond to 1350 generations ago (6_750 years ago).
    '''''
    for name, population in zip(vcf.samples, populations):
        if name == 'SAMN04028906':
            age = 1350
        else:
            age = 0
        samples.add_individual(ploidy=2, metadata={'ID':name}, population=population, time=age)

chromosome_length = {
    '1':158534110,
    '2':136231102,
    '3':121005158,
    '4':120000601,
    '5':120089316,
    '6':117806340,
    '7':110682743,
    '8':113319770,
    '9':105454467,
    '10':103308737,
    '11':106982474,
    '12':87216183,
    '13':83472345,
    '14':82403003,
    '15':85007780,
    '16':81013979,
    '17':73167244,
    '18':65820629,
    '19':63449741,
    '20':71974595,
    '21':69862954,
    '22':60773035,
    '23':52498615,
    '24':62317253,
    '25':42350435,
    '26':51992305,
    '27':45612108,
    '28':45940150,
    '29':51098607,
    'X':139009144,
    'Y':43300181,
    'M':16340
} 

  
#---------------------------------------------------------------------------------------------------------------------------  
# Inputs ----------------------------------------------------------
chromosome = str(sys.argv[1])
vcf_location = sys.argv[2]
metadata_location = sys.argv[3]
ancestral_location = sys.argv[4]

ne = sys.argv[5]

# As python objects ------------------------------------------------
vcf = cyvcf2.VCF(vcf_location, strict_gt=True)
name_vcf = vcf_location.split('/')[-1].strip('.vcf')

metadata = pd.read_csv(metadata_location, dtype=str)

ancestral_allele = pd.read_csv(ancestral_location, header=None, sep='\t')
#ancestral_allele = ancestral_allele[['pos', 'ances']]
#---------------------------------------------------------------------------------------------------------------------------

# Start inference
# (1) Generate sample files
#  creates the sampleData adding metadata
output = name_vcf + '.samples'
with tsinfer.SampleData(path=output, sequence_length=chromosome_length[str(chromosome)], num_flush_threads=10, max_file_size=2**30) as sample_file:
    populations = add_populations_bulls(vcf=vcf, samples=sample_file, metadata=metadata)
    add_diploid_individuals(vcf=vcf, samples=sample_file, populations=populations, metadata=metadata)
    add_diploid_sites(vcf=vcf, ancestral_state_file=ancestral_allele, sample_data=sample_file)

print(
    'Sample file  created for {} smaples'.format(sample_file.num_samples)
    + '({} individuals)'.format(sample_file.num_individuals)
    + 'with {} variable sites'.format(sample_file.num_sites)
)
# add age information (stored in metadata) to sampleData.individuals_time
#times = samples.individuals_time[:]
#copy = samples.copy(sample_file +  '.new', max_file_size=2**30)
#for indiv in samples.individuals():
#    times[indiv.id] = int(meta.Age[meta.ID == indiv.metadata['ID']])
#copy.individuals_time[:] = times
#copy.finalise()

#samples = copy

#samples = tsinfer.load(sample_file)

# (2) Generate ancestral samples and truncate eliminating long ancestrals
# REVIEWING VALUES FOR LOWER & UPPER CUTS
anc = tsinfer.generate_ancestors(
    sample_file,
    num_threads=10,
    progress_monitor=True,
).truncate_ancestors(
    lower_time_bound=0.4, upper_time_bound=0.6, length_multiplier=2,                 
)
# infer tree for ancestrals
inferred_anc_ts = tsinfer.match_ancestors(
    sample_file,
    anc,
    num_threads=10,
    recombination_rate=1.1e-8,
    mismatch_ratio=1,
    progress_monitor=True,
)
inferred_anc_ts.dump(name_vcf + '.atrees')

# (3) Infer final tree sequence by matching ancestrals tree and sampleData
inferred_ts = tsinfer.match_samples(
    sample_file,
    inferred_anc_ts,
    num_threads=10, 
    recombination_rate=1.1e-8,
    mismatch_ratio=1,
    progress_monitor=True,
).simplify(keep_unary=False)

inferred_ts.dump(name_vcf + '_final.trees')

# end of inferrence ----------------------------------------------------------------------

ne = 60_000
#ts_raw = tskit.load(f'{prefix}.trees')

ts_processed = tsdate.preprocess_ts(
    tree_sequence=inferred_ts,
    minimum_gap=1000000,
    remove_telomeres=True,
    filter_sites=False,
)
ts_processed.dump(f'{name_vcf}_final.processed.trees')

prior_obj = tsdate.build_prior_grid(
    tree_sequence=ts_processed, 
    timepoints=10,
    Ne=ne,
    approximate_priors=False,
    approx_prior_size=None,
    prior_distribution='lognorm',
    eps=5, # minimum distance separating time points in the time grid (specifies error factorn in time difference calculations)
    progress=True
)

# Process trees 
ts_dated = tsdate.date(
    ts_processed,
    time_units='generations',
    num_threads=10,
    mutation_rate=2.2e-9,
    priors=prior_obj,
    ignore_oldest_root=False,
    progress=True
)

ts_dated.dump(f'{name_vcf}_final.dated.trees')

