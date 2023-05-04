import itertools
import msprime
import numpy as np
import pandas as pd
import time
import allel
from scipy.sparse.linalg import LinearOperator, eigsh

from plotnine import *
%matplotlib inline

%reload_ext memory_profiler
%reload_ext line_profiler

def genetic_relatedness_matrix(ts, sample_sets, indexes, mode):
    n = len(sample_sets)
    K = np.zeros((n,n))
    K[np.triu_indices(n)] = ts.genetic_relatedness(sample_sets, indexes, mode = mode, proportion=False, span_normalise=False)
    K = K + np.triu(K,1).transpose()
    return K

def allelPCA(ts):
    haps = allel.HaplotypeArray(ts.genotype_matrix())
    gns = haps.to_genotypes(ploidy=2)
    allel_pc, fitted = allel.pca(gns.to_n_ref(), n_components=6, scaler=None)

    return allel_pc

def tsPCA(ts, mode = "branch", strategy = "direct"):

    # Set up samples
    n_ind = int(ts.num_samples / 2)
    sample_sets = [(2 * i, (2 * i) + 1) for i in range(n_ind)]
    W_samples = np.array([[float(u in A) for A in sample_sets] for u in ts.samples()])

    if strategy == "direct":
        # Direct genetic_relatedness matrix
        n = len(sample_sets)
        indexes = [
            (n1, n2) for n1, n2 in itertools.combinations_with_replacement(range(n), 2)
        ]
        K = genetic_relatedness_matrix(ts, sample_sets, indexes, mode)
        eigval_mat, eigvec_mat = eigsh(K)
        prin_comp = (K @ eigvec_mat[:, ::-1]) / np.sqrt(eigval_mat[::-1])
    else:
        # Linear operator
        indexes = [(i, n_ind) for i in range(len(sample_sets))]
        def mat_mul_stat(a):
            W = np.c_[W_samples, W_samples @ a]
            return ts.genetic_relatedness_weighted(
                W, indexes=indexes, mode=mode, span_normalise=False
                )
        A = LinearOperator((n_ind, n_ind), matvec=mat_mul_stat)
        eigval_linop, eigvec_linop = eigsh(A)
        prin_comp = (A @ eigvec_linop[:, ::-1]) / np.sqrt(eigval_linop[::-1])

    return prin_comp
