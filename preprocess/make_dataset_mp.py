import numpy as np
import os
import pandas as pd
import numpy as np
import argparse

import os.path
from multiprocessing import Process, Manager
import sys
import time
import itertools 
from tqdm import tqdm

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

cadd_anno = ["GC", "CpG", "SIFTcat", "SIFTval","PolyPhenCat",
             "PolyPhenVal", "bStatistic", "priPhCons","mamPhCons","verPhCons",
             "priPhyloP","mamPhyloP","verPhyloP","GerpN","GerpS","PHRED"]

# TODO are we interested in:
# 'non_coding_transcript_variant', 'frameshift_variant', or 'PolyPhenCat' 'Unknown' category?
cats = {
    # TODO these are strange, why does PolyPhen include 'benign' and 'unknown' as categories but SIFT doesn't
    "SIFTcat": np.array(["deleterious", "tolerated"]),
    "PolyPhenCat": np.array(["benign", "possibly_damaging", "unknown"]),
    "LoF": np.array(["LoF_HC", "LoF_LC"]),
    "Consequence": np.array(['3_prime_UTR_variant',
                     '5_prime_UTR_variant',
                     'TF_binding_site_variant',
                     'downstream_gene_variant',
                     'intergenic_variant',
                     'intron_variant',
                     'missense_variant',
#                      'non_coding_transcript_variant',
#                      'non_coding_transcript_exon_variant',
                     'regulatory_region_variant',
                     'splice_acceptor_variant',
                     'splice_donor_variant',
                     'splice_region_variant',
                     'stop_gained',
                     'synonymous_variant',
                     'upstream_gene_variant',
                     'frameshift_variant',
                    ]),
}


def to_cat(var, cat):
    if len(var) and cat == 'LoF': var = f'LoF_{var}'
    #if var == 'nan' and cat == 'PolyPhenCat': var = 'unknown'
    return pd.DataFrame(np.isin(cats[cat],var).astype(int).reshape(1,-1), columns=cats[cat])

def query_tabix(ch, s, e, SNV=True, indel=True):
    # Tabix indexing is 1-based, and the range (s-e) is inclusive of s,e.
    # VCF indexing is also 1-based
    tables = []
    if indel:
        # search starting at pos-1 since cadd files list deletions as {prev base}{ref} / {ref}
        os.system(f"tabix -h {cadd_file_indel} {ch}:{s-1}-{e} > /tmp/{ch}:{s}-{e}_indel.tsv")
        tables.append(pd.read_table(f"/tmp/{ch}:{s}-{e}_indel.tsv", header=1))
    if SNV: 
        os.system(f"tabix -h {cadd_file} {ch}:{s}-{e} > /tmp/{ch}:{s}-{e}.tsv")
        tables.append(pd.read_table(f"/tmp/{ch}:{s}-{e}.tsv", header=1))
    return pd.concat(tables)

def get_cadd(ch,pos,ref,alts):
    alt_to_cadd = {}
    ch = int(ch.strip('chr'))
    pos = int(pos)
    end = pos + len(ref) - 1

    alt_len = np.array([ (0 if alt == '*' else len(alt)) for alt in alts]) 
    SNV, indel = np.any(alt_len == 1), np.any(alt_len != 1)
    cadd_table = query_tabix(ch, pos, end, SNV, indel)

    if '*' in alts: # deletions are indexed at pos-1 in cadd files 
        del_cadd_table = cadd_table.loc[(cadd_table["Ref"].str.match(f".+{ref}$")) & (cadd_table["Pos"] == pos-1)]
    cadd_table = cadd_table.loc[(cadd_table["Ref"] == ref) & (cadd_table["Pos"] == pos)]

    for alt in alts:
        if alt == '*': 
            del_alts = del_cadd_table["Ref"].str.replace(f'{ref}$','')
            cadd_table_alt = del_cadd_table.loc[del_cadd_table["Alt"] == del_alts]
        else: 
            cadd_table_alt = cadd_table[cadd_table["Alt"] == alt][cadd_anno]

        if not len(cadd_table_alt):
            # a lot of indels not in gnomad, so suppress this warning for now 
            #print(f"did not find {ch} {pos} {ref} {alt} in {cadd_file} or {cadd_file_indel}")
            alt_to_cadd[alt] =  [""] * len(cadd_anno)
            continue
        
        cadd = cadd_table_alt.max(numeric_only=True).array.astype(str).tolist()
        svi = np.where(cadd_table_alt.columns == "SIFTval")[0][0]
        pvi = np.where(cadd_table_alt.columns == "PolyPhenVal")[0][0]

        # Replace SIFT with categoriacal
        sci = np.where(cadd_table_alt.columns == "SIFTcat")[0][0]
        sc = to_cat(cadd_table_alt["SIFTcat"].astype(str).iloc[0], "SIFTcat") # implicitly takes a maximum over the table
        sc = sc.values.reshape(-1).astype(str).tolist()
        cadd[sci] = sc[0]
        for e in range(len(sc[1:])): cadd.insert(sci+e,sc[e])

        # Replace PolyPhen with categorical    
        pci = np.where(cadd_table_alt.columns == "PolyPhenCat")[0][0] + len(sc[1:])
        pc = to_cat(cadd_table_alt["PolyPhenCat"].astype(str).iloc[0], "PolyPhenCat") # implicitly takes maximum over the table
        pc = pc.values.reshape(-1).astype(str).tolist()
        cadd[pci]=pc[0]
        for e in range(len(pc[1:])): cadd.insert(pci+e,pc[e])
        alt_to_cadd[alt] = cadd
    return alt_to_cadd

def get_vep(l):
    if 'CSQ=' not in l[7]: return None

    vep_infos = l[7].split('CSQ=')[1].split(',')
    vep_infos = [v.split('|') for v in vep_infos]

    vep_df = pd.DataFrame(vep_infos, columns=vep_fields)

    # convert to categorical
    vep_cat_df = pd.DataFrame()
    for v in vep_df.index:
        var = to_cat(vep_df.loc[v]["Consequence"], "Consequence")
        lof = to_cat(vep_df.loc[v]["LoF"], "LoF")
        cats = pd.concat([var, lof], axis=1)
        vep_cat_df = pd.concat([vep_cat_df, cats], axis=0)
    vep_cat_df.index = vep_df.loc[:,"Gene"]
 
    # aggregate multiple annotations over same gene using max
    for g in vep_cat_df.index.unique():
        if len(vep_cat_df.loc[[g]]) > 1:
            df2 = pd.DataFrame(vep_cat_df.loc[[g]].max(0)).T
            df2.index = [g]
            vep_cat_df = pd.concat([vep_cat_df.drop(g), df2])

    return vep_cat_df
    
def process_line(l, col_names):
    
    l = l.strip("\n").split()
    ldict = dict(zip(col_names,l))

    # Process INFO field
    linfo = {}
    missing = 0
    for e in ldict["INFO"].split(';'):
        kv = e.split('=')
        if len(kv) == 2: linfo[kv[0]] = kv[1]
        elif len(kv) == 0: continue 
        elif len(kv) == 1: 
            linfo[missing] = kv[0]
            missing += 1
        else: print(e); assert False
    lgtype = [e.split(":")[0].replace('|', '/') for e in l[9:]]
    col_ids = col_names[9:]    
    rare_idx = np.where(~np.isin(lgtype, ["0/0", "./."]))[0]
    rare_ids = [col_ids[i] for i in rare_idx]
 
    # check which alleles are present
    idx_to_alt = {}
    idx = 1
    alts = ldict["ALT"].split(",")
    alts_present = []
    for ac,alt in zip(linfo["AC"].split(","), alts):
        if ac == '0': 
            idx += 1
            continue
        idx_to_alt[str(idx)] = alt
        if  int(ac) > 0: alts_present.append(alt) 
        idx += 1

    vep_cat_df = get_vep(l)
    if vep_cat_df is None: return None

    ch,pos,ref = ldict["CHROM"], ldict["POS"], ldict["REF"]
    alt_to_cadd = get_cadd(ch,pos,ref,alts_present)
    
    out_lines = []

    # Create annotation line for each id, gene pair
    for ridx in set(rare_idx):
        rid = col_ids[ridx]  
        alleles = lgtype[ridx].split("/")
        for allele in np.unique(alleles):
            if allele == '0': continue
            alt = idx_to_alt[allele]
            for gene in vep_cat_df.index:
                out_line = [rid, gene, ch, pos]
                out_line.extend(vep_cat_df.loc[gene].array.astype(str).tolist())
                out_line.extend(alt_to_cadd[alt])
                out_lines.append("\t".join(out_line) + "\n")

    genes = vep_cat_df.index.array
    return "".join(out_lines), genes, rare_ids


def get_header():
    header = ["SubjectID", "GeneName", "Chromosome", "Position"]
    header.extend(cats["Consequence"])
    header.extend(cats["LoF"])

    cadd_anno_cols = []
    for c in cadd_anno:
        if c[-3:].lower() ==  "cat":
            cadd_anno_cols.extend([f"{c[:-3]}_{cc}" for cc in cats[c]])
        else: cadd_anno_cols.append(c)
    header.extend(cadd_anno_cols)
    return "\t".join(header) + "\n"


def do_work(in_queue, out_list, col_names):
    t1  = time.time()
    while True:
        if not in_queue.empty():
            item = in_queue.get(timeout=1)
            line_no, line = item
            # exit signal 
            if line == None:
                return
            result = process_line(line, col_names)
            out_list.append((line_no, result))
        # timeout
        elif time.time() - t1 > 1: return

if __name__ == "__main__":

    argParser = argparse.ArgumentParser()
    argParser.add_argument("--nw", default=32, type=int, help="number of workers (processes)")
    argParser.add_argument("--buf", default=100, type=int, help="buffer size per process")
    argParser.add_argument("--pop", default="AFR", type=str)
    argParser.add_argument("--postfix_in",
                            #default='hg38a.ID.ba.VEP.rare', 
                            default='30x.ID.VEP.bedtools.rare', 
                            type=str)
    argParser.add_argument("--postfix_out",
                            default='ws',
                            type=str)
    argParser.add_argument("--data_dir", default='/oak/stanford/groups/smontgom/erobb/data', type=str)

    args = argParser.parse_args()
    cadd_file = f"{args.data_dir}/cadd/whole_genome_SNVs_inclAnno.tsv.gz"
    cadd_file_indel = f"{args.data_dir}/cadd/gnomad.genomes.r3.0.indel_inclAnno.tsv.gz"

    #vcf_in = f"AF.all.{args.pop}.{args.postfix_in}.vcf"
    #tsv_out = f'AF.all.{args.pop}.{args.postfix_in}.{args.postfix_out}.tsv'
    vcf_in = f"all.{args.pop}.{args.postfix_in}.vcf"
    tsv_out = f'all.{args.pop}.{args.postfix_in}.{args.postfix_out}.tsv'
    gene_outliers = f'gene_outliers_{args.pop}.tsv'
    vcf_file = f"{args.data_dir}/vep/{vcf_in}"

    print(vcf_file)
    vcf = open(vcf_file, 'r')
    l = None
    # Process VCF header
    while True:
        prev_l = l
        hidx = vcf.tell()
        l = vcf.readline()
        if l[:14] == "##INFO=<ID=CSQ":
            vep_fields = l.split("Format: ")[1].strip('>"').split('|')
        if l[0] != '#':
            col_names = prev_l[1:-1].split()
            vcf.seek(hidx)
            break

    out = open(f"{args.data_dir}/watershed/{tsv_out}", 'w')
    # TODO move the gout creation to residuals.ipynb
    gout = open(f"{args.data_dir}/eoutliers/{gene_outliers}", 'w')

    header = get_header()
    out.write(header)

    # for debugging
    #while True:
    #    line = vcf.readline()
    #    result = process_line(line, col_names)
    #    print(result[0])
 
    os.system(f"wc -l {vcf_file} > /tmp/wc_{args.pop}")
    n = int(open(f"/tmp/wc_{args.pop}", "r").readlines()[0].split()[0])
    num_workers = args.nw
    lines_ps = args.buf

    manager = Manager()
    results = manager.list()
    work = manager.Queue(num_workers)

    pool = []
    for i in range(num_workers):
        p = Process(target=do_work, args=(work, results, col_names))
        p.start()
        pool.append(p)

    # produce data
    iters = itertools.chain(vcf, (None,)*num_workers)
    for num_and_line in tqdm(enumerate(iters), total=n):
        work.put(num_and_line)

        if num_and_line[0] % (num_workers * lines_ps) == 0 and num_and_line[0] != 0:
            # join processes
            for p in pool:
                p.join()

            # write the results
            for r in sorted(results):
                out_line, genes, rid = r[1]
                if out_line is not None and out_line != "":
                    out.write(out_line)
                    for gene in genes:
                        gout.write("\t".join([gene] + rid) + "\n")
            out.flush()
            gout.flush()

            # restart workers
            manager = Manager()
            results = manager.list()
            work = manager.Queue(num_workers)
            pool = []
            for i in range(num_workers):
                p = Process(target=do_work, args=(work, results, col_names))
                p.start()
                pool.append(p)

    out.close()  
    gout.close()
    print("finished processing vcf.")


