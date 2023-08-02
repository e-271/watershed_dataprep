import numpy as np
import os
import pandas as pd
import numpy as np
import argparse
import json

from collections import defaultdict


def sample_keys(d, n):
    kz =  np.random.choice([*d.keys()], size=n, replace=False)
    import pdb; pdb.set_trace()
    return {k: d[k] for k in kz}

def label_pairs(tsv_in, tsv_out, gene_outliers):
    agg_df = pd.read_table(tsv_in)
    agg_df = agg_df.set_index(["GeneName","SubjectID"], drop=False)
    pairs_all = json.load(open(gene_outliers, "r"))
    pairs_n2 = {k: pairs_all[k] for k in pairs_all if len(pairs_all[k]) == 2}
    pairs_ng2 = {k: pairs_all[k] for k in pairs_all if len(pairs_all[k]) > 2}

    # Divide out n >= 2 genes for "test set"
    num_test = int(len(pairs_all) * 0.3)
    if len(pairs_n2) < num_test:
        # May create a smaller test set if there are few N>2 genes
        pairs_test = dict(pairs_n2, **sample_keys(pairs_ng2, min(len(pairs_ng2), num_test-len(pairs_n2)))
    else:
        pairs_test = sample_keys(pairs_n2, num_test) 

    # replace z-scores with p-values
    def normalpdf(x): return np.exp(-x**2 / 2) / np.sqrt(2 * np.pi)
    agg_df["eOutlier"] = normalpdf(agg_df["eOutlier"])

    # drop af
    agg_df = agg_df.drop("AF", axis=1)
 
    # Replace null annotations with median
    for c in agg_df.columns[2:-1]:
        nmask = agg_df[c].isnull()
        # Drop columns that don't vary
        if agg_df[c].var() <  10e-10:
            agg_df = agg_df.drop(c, axis=1)
        elif nmask.any():
            med = agg_df.loc[np.invert(nmask), c].median()
            agg_df.loc[nmask, c] = med

    # Add N > 1 group IDs, and make a list of all n > 1 indices   
    agg_df["N2pair"] = np.nan
    pi = 0
    for k in pairs_test.keys():
        gene = k.split("_")[0]
        pairs = [(gene, sid) for sid in pairs_test[k]]
        pairs = agg_df.index.intersection(pairs)
        if len(pairs) < 2: continue
        pi += 1
        agg_df.loc[pairs[:2], "N2pair"] = pi
        if len(pairs) > 2:
            agg_df = agg_df.drop(pairs[2:])

    agg_df = agg_df.sort_values("N2pair") # assumed by WatershedR library!
    agg_df.loc[~agg_df["N2pair"].isnull(), "N2pair"] =  agg_df.loc[~agg_df["N2pair"].isnull(), "N2pair"].astype(int)
    agg_df.loc[agg_df["N2pair"].isnull(), "N2pair"] = "NA"

    agg_df.to_csv(tsv_out, index=False, sep="\t")



if __name__ == '__main__':
    argParser = argparse.ArgumentParser()
    argParser.add_argument("--pop", default="ESN", type=str)
    argParser.add_argument("--postfix_in",
                            #default='hg38a.ID.ba.VEP.rare.ws.gencode.phyloP.agg.eout', 
                            default='30x.ID.VEP.bedtools.rare.ws.gencode.phyloP.agg.eout',
                            type=str)
    argParser.add_argument("--postfix_out",
                            default='pairs',
                            type=str)
    argParser.add_argument("--data_dir", default='/oak/stanford/groups/smontgom/erobb/data', type=str)

    args = argParser.parse_args()
    pop = args.pop

    tsv_in = f'{args.data_dir}/watershed/all.{args.pop}.{args.postfix_in}.tsv'
    tsv_out = f'{args.data_dir}/watershed/all.{args.pop}.{args.postfix_in}.{args.postfix_out}.tsv'
    gene_outliers = f"{args.data_dir}/watershed/variant_pairs_{args.pop}.json"

    label_pairs(tsv_in, tsv_out, gene_outliers)



