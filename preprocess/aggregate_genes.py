import os
import pandas as pd
import numpy as np
import argparse
import warnings
from tqdm import tqdm

warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)

def aggregate(df):
    
    def _agg_func(x):
        if x.name in ['AF', 'distTSS', 'distTES']: return x.min()
        elif x.name == 'num_rare_variants': return x.sum()
        else: return x.max()

    def _apply_agg(df):
        df_agg = df.apply(_agg_func)
        return pd.DataFrame(df_agg).T
    
    k = 10e3
    cur_gene = None
    n = len(df)
    df = df.sort_values(["SubjectID", 'Chromosome', 'Position'])
    df = df.set_index(["SubjectID", "GeneName"], drop=False)
    df.insert(5, "num_rare_variants", 1)
    agg_df = pd.DataFrame(columns=df.columns)
    agg_df = agg_df.set_index(["SubjectID", "GeneName"], drop=False)
    i = 0
 
    pbar = tqdm(total=n) 
    while i < n:
        # Find start & end positions for current gene
        cur_id = df.iloc[i]["SubjectID"]
        cur_gene = df.iloc[i]["GeneName"]

        if (cur_id, cur_gene) in agg_df.index:
            i += 1
            pbar.update(1)
            continue

        gene_df = df.loc[(cur_id, cur_gene)]

        # Find min & max position for current gene
        ch = gene_df["Chromosome"].iloc[0] # should all be on same chromosome right?
        pos = gene_df["Position"]
        pmin, pmax = pos.min(), pos.max()

        # Find += 10kb windows
        astart, aend = i, i+1
        while astart > 0 and \
          df.iloc[astart-1]["Position"] >= pmin - k and \
          df.iloc[astart-1]["SubjectID"] == cur_id and \
          df.iloc[astart-1]["Chromosome"] == ch:
            astart -= 1
        while aend < n and \
          df.iloc[aend]["Position"] <= pmax + k and \
          df.iloc[aend]["SubjectID"] == cur_id and \
          df.iloc[aend]["Chromosome"] == ch:
            aend += 1
        assert aend >= i+1 and astart <= i

        # Aggregate over gene += 10kb
        agg = _apply_agg(df.iloc[astart:aend])
        agg["SubjectID"] = cur_id
        agg["GeneName"] = cur_gene
        import pdb; pdb.set_trace()
        agg = agg.set_index(["SubjectID", "GeneName"], drop=False)
        if (cur_id, cur_gene) in agg_df.index:
            print(i, cur_id, cur_gene)
            assert False
        agg_df = pd.concat([agg_df, agg])

        i += 1
        pbar.update(1)

    assert len(agg_df.index.unique()) == len(agg_df.index)
    agg_df = agg_df.drop(["Chromosome", "Position"], axis=1)

    return agg_df


if __name__ == '__main__':
    argParser = argparse.ArgumentParser()
    argParser.add_argument("--pop", default="ESN", type=str)
    argParser.add_argument("--data_dir", default='/oak/stanford/groups/smontgom/erobb/data', type=str)
    args = argParser.parse_args()

    tsv_in = f'AF.all.{args.pop}.hg38a.ID.ba.VEP.gencode.phyloP.rare.ws.tsv'
    tsv_out =  f'AF.all.{args.pop}.hg38a.ID.ba.VEP.gencode.phyloP.rare.agg.ws.tsv'
    tsv_file = f'{args.data_dir}/watershed/{tsv_in}'
    tsv_file_out = f'{args.data_dir}/watershed/{tsv_out}'    

    var_df = pd.read_table(tsv_file)
    agg_df = aggregate(var_df)
    #agg_df = agg_df.fillna("NA")
    agg_df.to_csv(tsv_file_out, sep="\t", index=False)
    

