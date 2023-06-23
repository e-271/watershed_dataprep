import os
import pandas as pd
import numpy as np
import argparse
import warnings
from tqdm import tqdm

warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)

def aggregate(df, genes_keep):
    
    def _agg_func(x):
        if x.name in ['AF', 'distTSS', 'distTES']: 
            return x.min()
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
    df = df.sort_index()
    agg_df = pd.DataFrame(columns=df.columns)
    agg_df = agg_df.set_index(["SubjectID", "GeneName"], drop=False)
    i = 0
 
    pbar = tqdm(total=n) 
    drop_genes = []
    while i < n:
        # Find start & end positions for current gene
        cur_id = df.iloc[i]["SubjectID"]
        cur_gene = df.iloc[i]["GeneName"]

        if (cur_id, cur_gene) in agg_df.index or cur_gene not in genes_keep: 
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


impute_values = {
"distTSS": np.nan, # These should disappear after aggregation I think? Valid NaN for variants not in a coding region?
"distTES": np.nan,
'AF': np.nan, 

'3_prime_UTR_variant': 0, 
'5_prime_UTR_variant': 0,
'TF_binding_site_variant': 0, 
'downstream_gene_variant': 0,
'intergenic_variant':0 , 
'intron_variant': 0, 
'missense_variant': 0,
'regulatory_region_variant': 0, 
'splice_acceptor_variant': 0,
'splice_donor_variant': 0, 
'splice_region_variant': 0, 
'stop_gained': 0,
'synonymous_variant': 0, 
'upstream_gene_variant': 0, 
'frameshift_variant': 0,
       
'LoF_HC': 0, 
'LoF_LC': 0, 

'GC': 0.418, 
'CpG': 0.024, 

'SIFT_deleterious': 0, 
'SIFT_tolerated': 0,
'SIFTval': 0, 

'PolyPhen_benign': 0, 
'PolyPhen_possibly_damaging': 0,      
'PolyPhen_unknown': 0, 
'PolyPhenVal': 0, 

'bStatistic': 800.261,
 
'priPhCons': 0,
'mamPhCons': 0, 
'verPhCons': 0,

"phyloP": 0,
'priPhyloP': -0.029, 
'mamPhyloP': -0.005, 
'verPhyloP': 0.042,
       
'GerpN': 3, 
'GerpS': -0.2, 
'PHRED': 0,

'eOutlier': np.nan
}

def impute_missing(df):
    for c in df.columns:
        if np.issubdtype(df[c].dtype, np.number): # numeric type column
            w = np.isnan(df[c])
            if not np.any(w): continue
            df.loc[w,c] = impute_values[c] 
        elif np.issubdtype(df[c].dtype, np.flexible): # string type column
            raise NotImplementedError
        elif np.issubdtype(df[c].dtype, np.bool_): # bool type column
            raise NotImplementedError
        elif np.issubdtype(df[c].dtype, np.object_): # object type
            w = pd.isnull(df[c])
            wstr = df[c].astype(str).str.fullmatch("(NA)|(NaN)|(NAN)")
            w = w | wstr
            if not np.any(w): continue
            df.loc[w,c] = impute_values[c]
        else: assert False
    return df


if __name__ == '__main__':
    argParser = argparse.ArgumentParser()
    argParser.add_argument("--pop", default="ESN", type=str)
    argParser.add_argument("--data_dir", default='/oak/stanford/groups/smontgom/erobb/data', type=str)
    args = argParser.parse_args()

    tsv_in = f'AF.all.{args.pop}.hg38a.ID.ba.VEP.gencode.phyloP.rare.ws.tsv'
    tsv_out =  f'AF.all.{args.pop}.hg38a.ID.ba.VEP.gencode.phyloP.rare.agg.ws.tsv'
    tsv_file = f'{args.data_dir}/watershed/{tsv_in}'
    tsv_file_out = f'{args.data_dir}/watershed/{tsv_out}'    
    gid_file = f'{args.data_dir}/gencode/gencode.v43.gene_ids.protein_lincRNA.txt'

    gids = pd.read_table(gid_file, header=None)
    gids = gids[0].str.split(".", expand=True)[0].values
    var_df = pd.read_table(tsv_file)
    var_df = impute_missing(var_df)

    agg_df = aggregate(var_df, gids)
    agg_df.to_csv(tsv_file_out, sep="\t", index=False)
    

