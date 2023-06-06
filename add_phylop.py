"""
Add phyloP annotations to Watershed dataset.

Prereqs:
  # Download phyloP BigWig file
  wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/phyloP100way/hg38.phyloP100way.bw
  # BigWig python tool
  pip install pyBigWig


"""

import pandas as pd
import numpy as np
import argparse
import os
import pyBigWig as pbw

def filelen(file):
    fname = file.split("/")[-1]
    os.system(f"wc -l {file} > /tmp/wc_{fname}")
    n = int(open(f"/tmp/wc_{fname}", "r").readlines()[0].split()[0])
    return n

def add_phyloP(df, phylop):
    df = df.sort_values(['Chromosome', 'Position'])
    df = df.set_index(["Chromosome", "Position"], drop=False)
    bw = pbw.open(phylop)
    
    def _query(row): 
        ch, pos = row.name
        return bw.values(f"chr{ch}", pos, pos+1)[0]

    pp = df.apply(_query, 
             axis=1, 
             result_type='expand')
    df.insert(4, "phyloP", pp)
    return df
    
if __name__ == '__main__':
    argParser = argparse.ArgumentParser()
    argParser.add_argument("--pop", default="ESN", type=str)
    argParser.add_argument("--data_dir", 
                            default='/oak/stanford/groups/smontgom/erobb/data/watershed', type=str)
    argParser.add_argument("--phylop", 
                            default='oak/stanford/groups/smontgom/erobb/data/hg38.phyloP100way.bw',
                            type=str)
    args = argParser.parse_args()

    tsv_in = f'AF.all.{args.pop}.hg38a.ID.ba.VEP.rare.ws.tsv'
    tsv_out =  f'AF.all.{args.pop}.hg38a.ID.ba.VEP.gencode.phyloP.rare.ws.tsv'
    tsv_file = f'{args.data_dir}/{tsv_in}'
    tsv_file_out = f'{args.data_dir}/{tsv_out}'

    var_df = pd.read_table(tsv_file)
    pp_df = add_phyloP(var_df, args.phylop)
    pp_df.to_csv(tsv_file_out, sep="\t", index=False)

