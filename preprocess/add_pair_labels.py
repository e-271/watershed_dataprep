import numpy as np
import os
import pandas as pd
import numpy as np
import argparse

from collections import defaultdict

argParser = argparse.ArgumentParser()
argParser.add_argument("--pop", default="AFR", type=str)
argParser.add_argument("--data_dir", default='/oak/stanford/groups/smontgom/erobb/data/watershed', type=str)

args = argParser.parse_args()
pop = args.pop
data_dir = args.data_dir

tsv_in = f'AF.all.{pop}.hg38a.ID.ba.VEP.rare.agg.ws.tsv'
tsv_train =  f'AF.all.{pop}.hg38a.ID.ba.VEP.rare.agg.ws.pairs.tsv'
print(f"{data_dir}/{tsv_train}")

agg_df = pd.read_table(f"{data_dir}/{tsv_in}")
agg_df = agg_df.set_index(["SubjectID","GeneName"], drop=False)

# Get list of outlier individuals per gene
gene_outliers = f'gene_outliers_{pop}.tsv'
gout = open(f"{data_dir}/{gene_outliers}", 'r')

# Make a dictionary of gene : outliers
genes = defaultdict(list)
for l in gout:
    l = l.strip("\n").split()
    genes[l[0]].extend(l[1:])

# Divide out n = 2 genes for "test set"
# TODO will want to change how the Watershed packages does this eventually i think.
genes_n2 = {k: genes[k] for k in genes if len(genes[k]) >= 2}

# TODO move this into make_dataset_mp.py
# Replace eOutlier binary score with original value & scale.
eout_file = f"/oak/stanford/groups/smontgom/erobb/data/watershed/{pop}_exprResiduals.tsv"
eout_keep_file = f"{data_dir}/ids_outlier_filtered_{pop}_t3f75.txt"
eout_df = pd.read_table(eout_file, sep="\t", index_col=0).T
# drop transcript information
gene_names = [c.split(".")[0] for c in eout_df.columns]
assert(len(np.unique(gene_names)) == len(eout_df.columns))
eout_df.columns = gene_names
inds_keep = pd.read_table(eout_keep_file, sep=" ", index_col=0, header=None).T
eout_df = eout_df.loc[eout_df.index.intersection(inds_keep.columns)]
for sid, gene in agg_df.index.values:
    agg_df.loc[(sid, gene), "eOutlier"] = eout_df.loc[sid, gene]

# replace z-scores with p-values
def normalpdf(x): return np.exp(-x**2 / 2) / np.sqrt(2 * np.pi)
agg_df["eOutlier"] = normalpdf(agg_df["eOutlier"])

#for sid in eout_df.index:
#    for gene in eout_df.columns:
#        if (sid, gene) in agg_df.index:
#            import pdb; pdb.set_trace()
#            agg_df.loc[(sid, gene), "eOutlier"] = eout_df.loc[sid, gene]

# Replace null annotations with median
for c in agg_df.columns[2:-1]:
    nmask = agg_df[c].isnull()
    # Drop columns that don't vary
    #if len(agg_df[c].unique()) == 1:
    if agg_df[c].var() <  10e-10:
        agg_df = agg_df.drop(c, 1)
    elif nmask.any():
        med = agg_df.loc[np.invert(nmask), c].median()
        agg_df.loc[nmask, c] = med

# Add N > 1 group IDs, and make a list of all n > 1 indices   
agg_df["N2pair"] = "NA"
pi = 0
for k in genes_n2.keys():
    pairs = [(sid, k) for sid in genes_n2[k]]
    pairs = agg_df.index.intersection(pairs)
    if len(pairs) >= 2: 
        pi += 1
        agg_df.loc[pairs[:2], "N2pair"] = f"{pi}"


# TODO need to report the true eOutlier value and pass in the threshold
#agg_df = agg_df.sort_values("N2pair")
#agg_df = agg_df.fillna("NA")


# TODO for testing purposes
#agg_df["eOutlier"] = np.random.random(size=88510)
agg_df.to_csv(f"{data_dir}/{tsv_train}", index=False, sep="\t")



