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


def correct_pvalues(pvalues):
    # TODO not totally sure I did this right, ask someone
    # Probability that the minimum of a set of clusters is less than our observed minimum, 
    # which is the probability of our null hypothesis
    mp = min(pvalues)
    pr = 1 - (1-mp)**len(pvalues)
    return pr 

def add_sOutliers(agg_df, sout_file, scluster_file, gencode_file):

    # Splicing cluster counts
    sout_df = pd.read_table(sout_file, sep="\t", index_col=0)
    # Splicing cluster definitions (intron:intron junctions)
    scl_df = pd.read_table(scluster_file, sep=":", skiprows=1, header=None, index_col=0, dtype={0: str, 1: int, 2:int, 3: str})
    scl_dict = {}
    for ch in ["X", "Y"] +  list(range(1,23)):
        scl_dict[f"chr{str(ch)}"] = scl_df.loc[str(ch)].sort_values(1) # sort by first junction position

    gencode_df = pd.read_table(gencode_file, sep="\t", header=None)
    gencode_df[0] = gencode_df[0].str.replace('[.].*', '', regex=True)
    gencode_df = gencode_df.set_index(0)

    agg_df["sOutlier"] = np.nan
    missing = set()
    nm = 0
    for sid, gene in tqdm(agg_df.index):
        # get clusters associated
        ch, p1, p2 = gencode_df.loc[gene]
        jidx =  np.where((scl_dict[ch][1] >= p1) & (scl_dict[ch][2] <= p2))
        if not len(jidx[0]): continue
        clusters = scl_dict[ch].iloc[jidx][3].unique()
        pvalues = []
        for c in clusters:
            try: pvalues.append(sout_df.loc[c, sid])
            except: 
                # some clusters 'miss' (print an error) or silently disappear in SPOT.
                # still figuring out why
                missing.add(c)
                continue 
        if not len(pvalues): continue
        agg_df.loc[(sid, gene), "sOutlier"] = correct_pvalues(pvalues)
        if len(missing) > nm: 
            nm = len(missing)
            print(nm)
    print(missing)
    return agg_df

if __name__ == '__main__':

    argParser = argparse.ArgumentParser()
    argParser.add_argument("--pop", default="ESN", type=str)
    argParser.add_argument("--postfix_in",
                            #default='hg38a.ID.ba.VEP.rare.ws.gencode.phyloP.agg', 
                            default='30x.ID.VEP.bedtools.rare.ws.gencode.phyloP.agg.eout',
                            type=str)
    argParser.add_argument("--postfix_out",
                            default='sout',
                            type=str)
    argParser.add_argument("--data_dir", default='/oak/stanford/groups/smontgom/erobb/data', type=str)
    args = argParser.parse_args()
        
    sout_file = f"{args.data_dir}/soutliers/afr_emperical_pvalue.txt"
    scluster_file = f"{args.data_dir}/soutliers/AFRIntrons.leafcutter.txt"
    tsv_in =  f'all.{args.pop}.{args.postfix_in}.tsv'
    tsv_out =  f'all.{args.pop}.{args.postfix_in}.{args.postfix_out}.tsv'
    tsv_file = f'{args.data_dir}/watershed/{tsv_in}'
    tsv_file_out = f'{args.data_dir}/watershed/{tsv_out}'
    gencode_file = f'{args.data_dir}/gencode/gencode.v43.gene_pos.tsv'

    agg_df = pd.read_table(tsv_file)
    agg_df = agg_df.set_index(["SubjectID", "GeneName"], drop=False)
    agg_df = agg_df.sort_index() # for speed

    df = add_sOutliers(agg_df, sout_file, scluster_file, gencode_file)
    df.to_csv(tsv_file_out, sep="\t", index=False)





