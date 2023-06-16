import os
import pandas as pd
import numpy as np
import argparse
from collections import defaultdict
from tqdm import tqdm

def add_pair_labels(df, genefile):
    gout = open(genefile, 'r')
    genes = defaultdict(list)
    for l in gout:
        l = l.strip("\n").split()
        genes[l[0]].extend(l[1:])
    genes_n2 = {k: genes[k] for k in genes if len(genes[k]) >= 2}
    df["N2pair"] = "NA"
    pi = 0
    for k in tqdm(genes_n2.keys()):
        pairs = [(sid, k) for sid in genes_n2[k]]
        pairs = df.index.intersection(pairs)
        if len(pairs) >= 2:
            pi += 1
            df.loc[pairs[:2], "N2pair"] = f"{pi}"
    return df
    

def replace_zscores_with_pvalues(df, col_name="eOutlier"):
    def normalpdf(x): return np.exp(-x**2 / 2) / np.sqrt(2 * np.pi)
    df[col_name] = normalpdf(df[col_name])
    return df
    
def replace_nulls(df):
    # Replace null annotations with median
    for c in df.columns[2:-1]:
        nmask = df[c].isnull()
        # Also drop columns with 0 variance since they will cause problems for Watershed
        if df[c].var() <  10e-10:
            df = df.drop(c, 1)
        elif nmask.any():
            med = df.loc[np.invert(nmask), c].median()
            df.loc[nmask, c] = med
    return df



if __name__ == '__main__':

    argParser = argparse.ArgumentParser()
    argParser.add_argument("--pop", default="ESN", type=str)
    argParser.add_argument("--data_dir", default='/oak/stanford/groups/smontgom/erobb/data/watershed', type=str)
    args = argParser.parse_args()        

    genefile = f'{args.data_dir}/gene_outliers_{args.pop}.tsv'
    tsv_in = f'AF.all.{args.pop}.hg38a.ID.ba.VEP.rare.agg.ws.tsv'
    tsv_train =  f'AF.all.{args.pop}.hg38a.ID.ba.VEP.rare.agg.ws.pairs.tsv'
    print(f"{args.data_dir}/{tsv_train}")

    agg_df = pd.read_table(f"{args.data_dir}/{tsv_in}")
    agg_df = agg_df.set_index(["SubjectID","GeneName"], drop=False)
    
    df = replace_nulls(agg_df)
    df = replace_zscores_with_pvalues(df)
    df = add_pair_labels(df, genefile)
    df.to_csv(f"{args.data_dir}/{tsv_train}", index=False, sep="\t")
    
    


