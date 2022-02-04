import tsdate
import tsinfer
import tskit

for chr in range(7,8):
    tree = tskit.load("Chr" + str(chr) + ".trees")

    newTree = tsdate.preprocess_ts(tree, remove_telomeres=True)

    newTree.dump("Chr" + str(chr) + "_processed.trees")