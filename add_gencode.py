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

def binsearch(df, pos, c1, c2):
    """
    Assumes dataframe is sorted by c1,c2.
    """
    n = len(df)
    if n == 0: return None
    if pos > df[c2].iloc[n//2]: # pos is to the right of "end"
        return binsearch(df.iloc[n//2+1:], pos, c1, c2)
    elif pos < df[c1].iloc[n//2]: # pos is to the left of "start"
        return binsearch(df.iloc[:n//2], pos, c1, c2)
    elif df[c1].iloc[n//2] <= pos and df[c2].iloc[n//2] >= pos:
        return df.iloc[n//2]
    else:
        return None

def add_gencode(df, gencode):
    gc = pd.read_table(gencode, header=None)
    gc = gc.sort_values([0,3,4])
    gc = gc.set_index([0],drop=False)
    
    def _query(row):
        ch, pos = row.index.unique()[0]
        gc_ch_pos = binsearch(gc.loc[f"chr{ch}"], pos, 3, 4)
        if gc_ch_pos is None: return [np.nan, np.nan]
        TSS = pos - gc_ch_pos[3] #gc.loc[w,3].values[0]
        TES = gc_ch_pos[4] - pos #gc.loc[w,4].values[0] - pos
        return [TSS, TES]

    ann = []
    for idx in tqdm(df.index):
        ann.append(_query(df.loc[idx]))
    ann = np.array(ann)

    df.insert(4, "distTSS", ann[:,0])
    df.insert(5, "distTES", ann[:,1])
    return df


if __name__ == '__main__':

    argParser = argparse.ArgumentParser()
    argParser.add_argument("--pop", default="ESN", type=str)
    argParser.add_argument("--data_dir", 
                            default='/oak/stanford/groups/smontgom/erobb/data/watershed',
                            type=str)
    argParser.add_argument("--gencode", 
                            default="/oak/stanford/groups/smontgom/erobb/data/watershed/gencode.v43.chr_patch_hapl_scaff.annotation.exons.protein_lincRNA.gtf", 
                            type=str)
    args = argParser.parse_args()
    
    tsv_in = f'AF.all.{args.pop}.hg38a.ID.ba.VEP.rare.ws.tsv'
    tsv_out =  f'AF.all.{args.pop}.hg38a.ID.ba.VEP.gencode.rare.ws.tsv'
    tsv_file = f'{args.data_dir}/{tsv_in}'
    tsv_file_out = f'{args.data_dir}/{tsv_out}'


    var_df = pd.read_table(tsv_file)
    var_df = var_df.sort_values(['Chromosome', 'Position'])
    var_df = var_df.set_index(["Chromosome", "Position"], drop=False)

    gdf = add_gencode(var_df, args.gencode)
    gdf.to_csv(tsv_file_out, sep="\t", index=False)





