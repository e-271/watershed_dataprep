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
        ch = ch.strip('chr')
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
                            default='/oak/stanford/groups/smontgom/erobb/data', type=str)
    argParser.add_argument("--postfix_in",
                            #default='hg38a.ID.ba.VEP.rare.ws.gencode', 
                            default='30x.ID.VEP.bedtools.rare.ws.gencode',
                            type=str)
    argParser.add_argument("--postfix_out",
                            default='phyloP',
                            type=str)
    argParser.add_argument("--phylop", 
                            default='hg38.phyloP100way.bw',
                            #default='241-mammalian-2020v2.bigWig',
                            type=str)
    args = argParser.parse_args()

    phylop_file = f'{args.data_dir}/phylop/{args.phylop}'
    tsv_in = f'all.{args.pop}.{args.postfix_in}.tsv'
    tsv_out =  f'all.{args.pop}.{args.postfix_in}.{args.postfix_out}.tsv'
    tsv_file = f'{args.data_dir}/watershed/{tsv_in}'
    tsv_file_out = f'{args.data_dir}/watershed/{tsv_out}'

    var_df = pd.read_table(tsv_file)
    pp_df = add_phyloP(var_df, phylop_file)
    pp_df.to_csv(tsv_file_out, sep="\t", index=False)

