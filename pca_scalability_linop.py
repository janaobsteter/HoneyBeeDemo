### Script to investigate scalability of various PCA methods

import tskit
import msprime
import itertools
import numpy as np
import pandas as pd
import allel
import time
from scipy.sparse.linalg import LinearOperator, eigsh


chromosome = 16
ts = tskit.load("Chr" + str(chromosome) + ".trees")
n_ind = 691
def simul(ts, seed, n_ind, mutation_rate):
    ts = msprime.simulate(
        n_ind * 2, 
        length=1e6, 
        recombination_rate=1e-8, 
        Ne=1e4, 
        mutation_rate=mutation_rate, 
        random_seed=seed
        )
    #n_ind = len(ts.samples())
    sample_sets = [(2 * i, (2 * i) + 1) for i in range(n_ind)]
    W_samples = np.array([[float(u in A) for A in sample_sets] for u in ts.samples()])
    indexes = [(i, n_ind) for i in range(len(sample_sets))]
    
    def mat_mul_stat(a):
        W = np.c_[W_samples, W_samples @ a]
        return ts.genetic_relatedness_weighted(
            W, indexes=indexes, mode="site", span_normalise=False
            )

    # Linear operator
    start_time = time.time()
    A = LinearOperator((n_ind, n_ind), matvec=mat_mul_stat)
    eigval_linop, eigvec_linop = eigsh(A)
    linop_pc = (A @ eigvec_linop[:, ::-1]) / np.sqrt(eigval_linop[::-1])
    end_time = time.time()
    time_linop = end_time - start_time

    # Direct genotype matrix
 
    haps = allel.HaplotypeArray(ts.genotype_matrix())
    gns = haps.to_genotypes(ploidy=2).to_n_ref()
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
         'mut_rate':[mutation_rate, mutation_rate],
         'n_ind':[n_ind, n_ind],
         'seed':[seed, seed],
         'p':[ts.num_sites, ts.num_sites],
         'time':[time_allel, time_linop],
         'diff':[0, diff_linop]}
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