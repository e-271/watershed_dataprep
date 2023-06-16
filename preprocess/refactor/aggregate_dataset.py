import numpy as np
import os
import pandas as pd
import numpy as np
import argparse


argParser = argparse.ArgumentParser()
argParser.add_argument("--pop", default="AFR", type=str)
argParser.add_argument("--data_dir", default='/oak/stanford/groups/smontgom/erobb/data/watershed', type=str)

args = argParser.parse_args()
pop = args.pop
data_dir = args.data_dir

data_dir = '/oak/stanford/groups/smontgom/erobb/data/watershed'
tsv_in = f'AF.all.{pop}.hg38a.ID.ba.VEP.rare.ws.tsv'
tsv_out =  f'AF.all.{pop}.hg38a.ID.ba.VEP.rare.agg.ws.tsv'
tsv_file = f'{data_dir}/{tsv_in}'
tsv_file_out = f'{data_dir}/{tsv_out}'


import warnings
warnings.filterwarnings("ignore")

def _agg(x):
    if x.name == 'af': return x.min()
    elif x.name == 'num_rare_variants': return x.sum()
    # TODO we might be interested to see the sum for other things, 
    # eg # missense variants in a gene?
    else: return x.max()

def aggregate(df):
    df_agg = df.apply(_agg)
    return pd.DataFrame(df_agg).T

k = 10e3

var_df = pd.read_table(tsv_file)
var_df = var_df.sort_values(["SubjectID", 'Chromosome', 'Position'])
var_df = var_df.set_index(["SubjectID", "GeneName"], drop=False)

cur_gene = None
n = len(var_df)
agg_df = pd.DataFrame(columns=var_df.columns)
agg_df = agg_df.set_index(["SubjectID", "GeneName"], drop=False)
i = 0
print(tsv_file_out )
while i < n:

    # Find start & end positions for current gene
    cur_id = var_df.iloc[i]["SubjectID"]
    cur_gene = var_df.iloc[i]["GeneName"]
        
    if (cur_id, cur_gene) in agg_df.index:
        i += 1 
        continue 

    gene_df = var_df.loc[(cur_id, cur_gene)]
    
    # Find min & max position for current gene
    ch = gene_df["Chromosome"].iloc[0] # should all be on same chromosome right?
    pos = gene_df["Position"]
    pmin, pmax = pos.min(), pos.max()
    
    # Find += 10kb windows
    astart, aend = i, i+1
    while astart > 0 and \
      var_df.iloc[astart-1]["Position"] >= pmin - k and \
      var_df.iloc[astart-1]["SubjectID"] == cur_id and \
      var_df.iloc[astart-1]["Chromosome"] == ch: 
        astart -= 1
    while aend < n and \
      var_df.iloc[aend]["Position"] <= pmax + k and \
      var_df.iloc[aend]["SubjectID"] == cur_id and \
      var_df.iloc[aend]["Chromosome"] == ch: 
        aend += 1
    assert aend >= i+1 and astart <= i
    
    # Aggregate over gene += 10kb
    agg = aggregate(var_df.iloc[astart:aend])
    agg["SubjectID"] = cur_id
    agg["GeneName"] = cur_gene
    agg = agg.set_index(["SubjectID", "GeneName"], drop=False)
    if (cur_id, cur_gene) in agg_df.index:
        print(i, cur_id, cur_gene)
        assert False
    agg_df = pd.concat([agg_df, agg])
    
    i += 1 
    if i % 1000 == 0: print(f"{i} / {n}")

assert len(agg_df.index.unique()) == len(agg_df.index)
agg_df = agg_df.drop(["Chromosome", "Position"], axis=1)
agg_df = agg_df.fillna("NA")
agg_df.to_csv(tsv_file_out, sep="\t", index=False)








