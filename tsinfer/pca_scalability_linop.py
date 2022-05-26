### Script to investigate scalability of various PCA methods

#import tskit
import itertools
import numpy as np
import pandas as pd
import allel
import time
import sys
from scipy.sparse.linalg import LinearOperator, eigsh
#sys.path.append("/home/x/JANA/github/PeterRalpha_tskit/python/")
import tskit

def genetic_relatedness_matrix(ts, sample_sets, mode):
    n = len(sample_sets)
    indexes = [(n1, n2) for n1, n2 in itertools.combinations_with_replacement(range(n), 2)]
    K = np.zeros((n,n))
    K[np.triu_indices(n)] = ts.genetic_relatedness(sample_sets, indexes, mode = mode, proportion=False, span_normalise=False)
    K = K + np.triu(K,1).transpose()
    return K


def simul(seed, n_ind, mutation_rate):

    def mat_mul_stat(a):
        W = np.c_[W_samples, W_samples @ a]
        return ts.genetic_relatedness_weighted(
            W, indexes=indexes, mode="site", span_normalise=False
            )

    ts = tskit.load("/media/x/0cfdc498-0e95-4f62-87a9-4c0c9f10b126/jana/Documents/1Projects/HoneybeeDemography/tsinfer/HiFiTrees/Chr3.trees")
    n_ind = ts.num_individuals
    sample_sets = [2 * i for i in range(n_ind)]
    W_samples = np.array([[float(u == A) for A in sample_sets] for u in ts.samples()]) #float in A for diploid
    indexes = [(i, n_ind) for i in range(len(sample_sets))]

    # Linear operator
    start_time = time.time()
    A = LinearOperator((n_ind, n_ind), matvec=mat_mul_stat)
    eigval_linop, eigvec_linop = eigsh(A)
    linop_pc = (A @ eigvec_linop[:, ::-1]) / np.sqrt(eigval_linop[::-1])
    end_time = time.time()
    time_linop = end_time - start_time

    # Direct genotype matrix
 
    haps = allel.HaplotypeArray(ts.genotype_matrix())
    gns = haps.to_genotypes(ploidy=1).to_n_ref()
    start_time = time.time()
    allel_pc = allel.pca(gns, n_components=6, scaler=None)[0]
    end_time = time.time()
    time_allel = end_time - start_time

    # Assess output
    diff_linop = np.max(
        np.abs(allel_pc) - np.abs(linop_pc)
    )

    out_df = pd.DataFrame(
        {'method':['allel', 'linop'],
         #'mut_rate':[mutation_rate, mutation_rate],
         'n_ind':[n_ind, n_ind],
         #'seed':[seed, seed],
         'p':[ts.num_sites, ts.num_sites],
         #'time':[time_allel, time_linop],
         #'diff':[0, diff_linop]
        }
    )
    return out_df

sim_df = pd.DataFrame({
    'method':[], 'mut_rate':[], 'n_ind':[], 'seed':[],
    'p':[], 'time':[], 'diff':[]
})
for seed in range(1, 11):
    for mutation_rate in [1e-6, 1e-7, 1e-8, 1e-9]:
        for n_ind in [8, 16, 32, 64, 128, 256, 512, 1024, 2048]:
            sim = simul(seed, n_ind, mutation_rate)
            sim_df = sim_df.append(sim, sort = True)

print(sim_df.to_string())
sim_df.to_csv("simul_linop.csv")











##################################
import json
import pandas as pd

sample_nodes = [ts.node(n) for n in ts.samples()]
sample_names = [
    json.loads(ts.individual(n.individual).metadata)["name"]
    for n in sample_nodes
]
pd.DataFrame({"Name": sample_names}).to_csv("SampleNames.csv")