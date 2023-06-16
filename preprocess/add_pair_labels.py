import numpy as np
import os
import pandas as pd
import numpy as np
import argparse

from collections import defaultdict


def label_pairs(tsv_in, tsv_out, gene_outliers):
    agg_df = pd.read_table(tsv_in)
    agg_df = agg_df.set_index(["SubjectID","GeneName"], drop=False)
    gout = open(gene_outliers)

    # Make a dictionary of gene : outliers
    genes = defaultdict(list)
    for l in gout:
        l = l.strip("\n").split()
        genes[l[0]].extend(l[1:])

    # Divide out n = 2 genes for "test set"
    genes_n2 = {k: genes[k] for k in genes if len(genes[k]) == 2}

    # replace z-scores with p-values
    def normalpdf(x): return np.exp(-x**2 / 2) / np.sqrt(2 * np.pi)
    agg_df["eOutlier"] = normalpdf(agg_df["eOutlier"])

    # Replace null annotations with median
    for c in agg_df.columns[2:-1]:
        nmask = agg_df[c].isnull()
        # Drop columns that don't vary
        if agg_df[c].var() <  10e-10:
            agg_df = agg_df.drop(c, 1)
        elif nmask.any():
            med = agg_df.loc[np.invert(nmask), c].median()
            agg_df.loc[nmask, c] = med

    # Add N > 1 group IDs, and make a list of all n > 1 indices   
    # TODO Taibo says to remove any additional entries for genes with n>2
    agg_df["N2pair"] = "NA"
    pi = 0
    import pdb; pdb.set_trace()
    for k in genes_n2.keys():
        pairs = [(sid, k) for sid in genes_n2[k]]
        pairs = agg_df.index.intersection(pairs)
        if len(pairs) >= 2: 
            pi += 1
            agg_df.loc[pairs[:2], "N2pair"] = f"{pi}"


    agg_df.to_csv(tsv_out, index=False, sep="\t")



if __name__ == '__main__':
    argParser = argparse.ArgumentParser()
    argParser.add_argument("--pop", default="ESN", type=str)
    argParser.add_argument("--data_dir", default='/oak/stanford/groups/smontgom/erobb/data', type=str)

    args = argParser.parse_args()
    pop = args.pop

    tsv_in = f'{args.data_dir}/watershed/AF.all.{args.pop}.hg38a.ID.ba.VEP.gencode.phyloP.rare.agg.ws.tsv'
    tsv_out = f'{args.data_dir}/watershed/AF.all.{args.pop}.hg38a.ID.ba.VEP.gencode.phyloP.rare.agg.ws.pairs.tsv'
    gene_outliers = f"{args.data_dir}/eoutliers/gene_outliers_{args.pop}.tsv"

    label_pairs(tsv_in, tsv_out, gene_outliers)



