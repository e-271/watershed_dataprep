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

def add_gencode(df, gencode):
    gc = pd.read_table(gencode, header=None)
    ann_df = gc[8].str.split(";", expand=True)
    gc["gene_id"] = ann_df[0].str.strip("gene_id ").str.strip('"')
    gc["gene_id"] = gc["gene_id"].str.split(".", expand=True)[0] 
    gc["transcript_id"]  = ann_df[1].str.strip("transcript_id ").str.strip('"')
    gc = gc.set_index(["gene_id"],drop=False)
    gc = gc.sort_index()    

    ann = []
    drop = []
    for idx in tqdm(df.index):
        row =  df.loc[idx]
        genes = row["GeneName"].unique()
        TSS, TES = [], []
        for gene in genes:
            if gene in gc.index:
                gc.loc[gene]
                ch, pos = row.index.unique()[0]
                tss, tes = gc.loc[gene, 3].min(), gc.loc[gene, 4].max()
                dTSS = abs(pos - tss)
                dTES = abs(tes - pos)
                TSS.append(dTSS)
                TES.append(dTES)
            else:
                TSS.append(np.nan)
                TES.append(np.nan)
        dTSS = min(TSS)
        dTES = min(TES)
        ann.append([dTSS, dTES])

    ann = np.array(ann)

    df.insert(4, "distTSS", ann[:,0])
    df.insert(5, "distTES", ann[:,1])
    return df


if __name__ == '__main__':

    argParser = argparse.ArgumentParser()
    argParser.add_argument("--pop", default="ESN", type=str)
    argParser.add_argument("--data_dir", default='/oak/stanford/groups/smontgom/erobb/data', type=str)
    args = argParser.parse_args()
        
    gencode_file = f'{args.data_dir}/gencode/gencode.v43.chr_patch_hapl_scaff.annotation.transcripts.protein_lincRNA.gtf'
    tsv_in = f'AF.all.{args.pop}.hg38a.ID.ba.VEP.rare.ws.tsv'
    tsv_out =  f'AF.all.{args.pop}.hg38a.ID.ba.VEP.gencode.rare.ws.tsv'
    tsv_file = f'{args.data_dir}/watershed/{tsv_in}'
    tsv_file_out = f'{args.data_dir}/watershed/{tsv_out}'


    var_df = pd.read_table(tsv_file)
    #var_df = var_df.sort_values(['Chromosome', 'Position'])
    var_df = var_df.set_index(["Chromosome", "Position"], drop=False)
    var_df = var_df.sort_index() # for speed

    gdf = add_gencode(var_df, gencode_file)
    gdf.to_csv(tsv_file_out, sep="\t", index=False)





