"""
Add gencode annotations to Watershed dataset.

Prepreqs:
  wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/gencode.v43.annotation.gtf.gz
  gunzip gencode.v43.annotation.gtf.gz
  sh scripts/process_gencode.sh

Output columns
  distTSS: distance to transcription start site
  distTES: distance to transcription end site
"""
import pandas as pd
import numpy as np
import argparse
import os
from tqdm import tqdm


def add_eOutliers(agg_df, eout_file, eout_keep_file):

    # Expression outlier scores from residuals file
    eout_df = pd.read_table(eout_file, sep="\t", index_col=0).T
    # drop transcript information
    gene_names = [c.split(".")[0] for c in eout_df.columns]
    assert(len(np.unique(gene_names)) == len(eout_df.columns))
    eout_df.columns = gene_names
    inds_keep = pd.read_table(eout_keep_file, sep=" ", index_col=0, header=None).T
    eout_df = eout_df.loc[eout_df.index.intersection(inds_keep.columns)]

    agg_df["eOutlier"] = np.nan
    for sid, gene in tqdm(agg_df.index):
        if gene in eout_df and sid in eout_df.index: 
            agg_df.loc[(sid,gene), "eOutlier"] = eout_df.loc[sid,gene]
    agg_df = agg_df[~agg_df["eOutlier"].isnull()]
    return agg_df

if __name__ == '__main__':

    argParser = argparse.ArgumentParser()
    argParser.add_argument("--pop", default="ESN", type=str)
    argParser.add_argument("--postfix_in",
                            #default='hg38a.ID.ba.VEP.rare.ws.gencode.phyloP.agg', 
                            default='30x.ID.VEP.bedtools.rare.ws.gencode.phyloP.agg',
                            type=str)
    argParser.add_argument("--postfix_out",
                            default='eout',
                            type=str)
    argParser.add_argument("--data_dir", default='/oak/stanford/groups/smontgom/erobb/data', type=str)
    args = argParser.parse_args()
        
    eout_file = f"{args.data_dir}/eoutliers/{args.pop}_exprResiduals.tsv"
    eout_keep_file = f"{args.data_dir}/eoutliers/ids_outlier_filtered_{args.pop}_t3.00f75.txt"
    tsv_in =  f'all.{args.pop}.{args.postfix_in}.tsv'
    tsv_out =  f'all.{args.pop}.{args.postfix_in}.{args.postfix_out}.tsv'
    tsv_file = f'{args.data_dir}/watershed/{tsv_in}'
    tsv_file_out = f'{args.data_dir}/watershed/{tsv_out}'

    agg_df = pd.read_table(tsv_file)
    agg_df = agg_df.set_index(["SubjectID", "GeneName"], drop=False)
    agg_df = agg_df.sort_index() # for speed

    df = add_eOutliers(agg_df, eout_file, eout_keep_file)
    df.to_csv(tsv_file_out, sep="\t", index=False)





