import json
import itertools
import numpy as np
import pandas as pd
import tskit
import time
#import allel
from scipy.sparse.linalg import LinearOperator, eigsh
import sys

chromosome = sys.argv[1]

print("Defining functions")
def genetic_relatedness_matrix(ts, sample_sets, indexes, mode):
    n = len(sample_sets)
    K = np.zeros((n, n))
    K[np.triu_indices(n)] = ts.genetic_relatedness(sample_sets, indexes, mode=mode, proportion=False,
                                                   span_normalise=False)
    K = K + np.triu(K, 1).transpose()
    return K


def allelPCA(ts):
    haps = allel.HaplotypeArray(ts.genotype_matrix())
    gns = haps.to_genotypes(ploidy=1)
    allel_pc, fitted = allel.pca(gns.to_n_ref(), n_components=6, scaler=None)

    return allel_pc


def tsPCA(ts, mode="branch", strategy="direct"):
    # Set up samples
    n_ind = int(ts.num_samples)
    sample_sets = [[i] for i in range(n_ind)]
    W_samples = np.array([[float(u == A) for A in sample_sets] for u in ts.samples()])

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

print("Loading the tree")
ts = tskit.load("Chr" + str(chromosome) + "_processed.trees")
print("Computing PCA")
print("Allele")
pc_allel = allelPCA(ts)
print("Site direct")
pc_site_direct = tsPCA(ts, 'site', 'direct')
print("Branch direct")
pc_branch_direct = tsPCA(ts,  'branch', 'direct')
pc_allel.shape
pc_site_direct.shape
pc_branch_direct.shape
print("Site linop")
pc_site_linop = tsPCA(ts,  'site', 'linop')
print("Branch linop")
pc_branch_linop = tsPCA(ts,  'branch', 'linop')

print("Getting metadata")
# Get metadata
sample_nodes = [ts.node(n) for n in ts.samples()]
sample_subspecies = [
    json.loads(ts.population(n.population).metadata)["subspecies"]
    for n in sample_nodes
]

# Get sample population
sample_pops = [
    json.loads(ts.population(n.population).metadata)["lineage"]
    for n in sample_nodes
]

sample_id_subspecies = [
    (n.id, json.loads(ts.population(n.population).metadata)["subspecies"])
    for n in sample_nodes
]

nonHNodes = [id for (id, subsp) in sample_id_subspecies if subsp != "hybrid"]

print("Writing pandas")
dfAllel = pd.DataFrame(pc_allel)
dfAllel.loc[:, "Subspecies"] = sample_subspecies
dfAllel.loc[:, "Lineage"] = sample_pops
dfAllel.to_csv("PCAAllele.csv")

dfSiteDirect = pd.DataFrame(pc_site_direct)
dfSiteDirect.loc[:, "Subspecies"] = sample_subspecies
dfSiteDirect.loc[:, "Lineage"] = sample_pops
dfSiteDirect.to_csv("PCASiteDirect.csv")

dfBranchDirect = pd.DataFrame(pc_branch_direct)
dfBranchDirect.loc[:, "Subspecies"] = sample_subspecies
dfBranchDirect.loc[:, "Lineage"] = sample_pops
dfBranchDirect.to_csv("PCABranchDirect.csv")

dfSiteLinop = pd.DataFrame(pc_site_linop)
dfSiteLinop.loc[:, "Subspecies"] = sample_subspecies
dfSiteLinop.loc[:, "Lineage"] = sample_pops
dfSiteLinop.to_csv("PCASiteLinop.csv")

dfBranchLinop = pd.DataFrame(pc_branch_linop)
dfBranchLinop.loc[:, "Subspecies"] = sample_subspecies
dfBranchLinop.loc[:, "Lineage"] = sample_pops
dfBranchLinop.to_csv("PCABranchLinop.csv")


ts = ts.simplify(nonHNodes)

colours = {
        "capensis": "#f7e225",
        "scutellata": "#feb72d",
        "unicolor": "#f07f4f",
        "ligustica": "#310597",
        "carnica": "#ad2793",
        "mellifera": "#228c8d",
        "caucasica": "#a8db34",
        "hybrid": "#83888f"
    }

coloursL = {
        "A": "#f7e225",
        "C": "#ad2793",
        "M": "#228c8d",
        "O": "#a8db34",
        "H": "#FFFFFF"
    }


colours_for_node = {}
for n in ts.samples():
    population_data = ts.population(ts.node(n).population)
    colours_for_node[n] = colours[json.loads(population_data.metadata)["subspecies"]]

colours_for_node_lineage = {}
for n in ts.samples():
    population_data = ts.population(ts.node(n).population)
    colours_for_node_lineage[n] = coloursL[json.loads(population_data.metadata)["lineage"]]


individual_for_node = {}
for n in ts.samples():
    individual_data = ts.individual(ts.node(n).individual)
    individual_for_node[n] = json.loads(individual_data.metadata)["name"]

### Csd tree
ts3 = tskit.load("/home/janao/Documents/1Projects/HoneybeeDemography/tsinfer/Chr3_processed.trees")
ts3 = ts3.simplify(nonHNodes)
csdTree = ts3.at(11775000)


csdTree.draw(
    path="csdTree.svg",
    height=700,
    width=1200,
    node_labels={},
    node_colours=colours_for_node,
)


# Hox gene tree
ts16 = tskit.load("/home/janao/Documents/1Projects/HoneybeeDemography/tsinfer/Chr16_processed.trees")
hoxtree = ts16.at(3885000)
hoxtree.draw(
    path="hoxTree.svg",
    height=700,
    width=1200,
    node_labels=individual_for_node,
    node_colours=colours_for_node,
)

# Yellowy gene tree
ts10 = tskit.load("/home/janao/Documents/1Projects/HoneybeeDemography/tsinfer/HiFiTrees/Chr10.trees")
ts10 = ts10.simplify(nonHNodes)
yellowtree = ts10.at(11506337)
yellowtree.draw(
    path="yellowtreeStart.svg",
    height=700,
    width=1200,
    node_labels={},
    node_colours=colours_for_node,
)


yellowtreeEnd = ts10.at(11513040)
yellowtree.draw(
    path="yellowtreeEnd.svg",
    height=700,
    width=1200,
    node_labels=individual_for_node,
    node_colours=colours_for_node,
)


yellowtree = ts10.at(11506337)
yellowtree.draw(
    path="yellowtreeStart_L.svg",
    height=700,
    width=1200,
    node_labels=individual_for_node,
    node_colours=colours_for_node_lineage,
)

yellowtreeEnd = ts10.at(11513040)
yellowtree.draw(
    path="yellowtreeEnd_L.svg",
    height=700,
    width=1200,
    node_labels=individual_for_node,
    node_colours=colours_for_node_lineage,
)

# PCA on the one tree
print("Site direct")
pc_site_direct_Csd = tsPCA(ts3, 'site', 'direct')
pc_site_direct_Yellow = tsPCA(ts10, 'site', 'direct')


sample_nodes = [ts10.node(n) for n in ts10.samples()]
sample_subspecies = [
    json.loads(ts10.population(n.population).metadata)["subspecies"]
    for n in sample_nodes
]

# Get sample population
sample_pops = [
    json.loads(ts10.population(n.population).metadata)["lineage"]
    for n in sample_nodes
]

dfSiteDirect = pd.DataFrame(pc_site_direct_Yellow)
dfSiteDirect.loc[:, "Subspecies"] = sample_subspecies
dfSiteDirect.loc[:, "Lineage"] = sample_pops
dfSiteDirect.to_csv("PCASiteDirect_Yellow.csv")

dfSiteDirect = pd.DataFrame(pc_site_direct_Csd)
dfSiteDirect.loc[:, "Subspecies"] = sample_subspecies
dfSiteDirect.loc[:, "Lineage"] = sample_pops
dfSiteDirect.to_csv("PCASiteDirect_Csd.csv")


# Plot with drw_SVG
styles = []
# Create a style for each population, programmatically (or just type the string by hand)
for colour, p in zip(["#ad2793", "#228c8d", "#feb72d", "#a8db34", "#f07f4f", "#f7e225","#310597"], ts10.populations()):
    # target the symbols only (class "sym")
    s = f".node.p{p.id} > .sym " + "{" + f"fill: {colour}" + "}"
    styles.append(s)
    print(f'"{s}" applies to nodes from population {json.loads(p.metadata)["subspecies"]} (id {p.id})')
hide_internal_symlabs = ".node:not(.leaf) > .sym, .node:not(.leaf) > .lab {visibility: hidden}"
css_string = " ".join(styles)
css_string = css_string + hide_internal_symlabs
print(f'CSS string applied:\n    "{css_string}"')


yellowtree.draw_svg(
    path="PCAYellow_SVG.svg",
    size=(700, 400),
    node_labels={},    # Remove all node labels for a clearer viz
    style=css_string,
    symbol_size=10 # Apply the stylesheet
)

csdTree.draw_svg(
    path="PCACsd_SVG.svg",
    size=(700, 400),
    node_labels={},    # Remove all node labels for a clearer viz
    style=css_string,
    symbol_size=10 # Apply the stylesheet
)

