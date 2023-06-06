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


argParser = argparse.ArgumentParser()
argParser.add_argument("--nw", default=32, type=int, help="number of workers (processes)")
argParser.add_argument("--buf", default=100, type=int, help="buffer size per process")
argParser.add_argument("--pop", default="AFR", type=str)
argParser.add_argument("--vep_dir", default='/oak/stanford/groups/smontgom/erobb/data/vep', type=str)
argParser.add_argument("--data_dir", default='/oak/stanford/groups/smontgom/erobb/data/watershed', type=str)
argParser.add_argument("--cadd", default="/oak/stanford/groups/smontgom/erobb/data/watershed/whole_genome_SNVs_inclAnno.tsv.gz", type=str)
argParser.add_argument("--gencode", default="/oak/stanford/groups/smontgom/erobb/data/watershed/gencode.v43.chr_patch_hapl_scaff.annotation.exons.protein_lincRNA.gtf", type=str)

args = argParser.parse_args()
pop = args.pop
vep_dir, data_dir = args.vep_dir, args.data_dir
CADD_FILE = args.cadd


vcf_in = f"AF.all.{pop}.hg38a.ID.ba.VEP.rare.vcf"
tsv_out = f'AF.all.{pop}.hg38a.ID.ba.VEP.rare.ws.tsv'
gene_outliers = f'gene_outliers_{pop}.tsv'
vcf_file = f"{vep_dir}/{vcf_in}"
#eout_file = f"{data_dir}/eOutlier_scores_{pop}_t3.txt"
eout_file = f"/oak/stanford/groups/smontgom/erobb/data/watershed/{pop}_exprResiduals.tsv"
eout_keep_file = f"{data_dir}/ids_outlier_filtered_{pop}_t3f75.txt"


vep_fields = 'Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|SYMBOL_SOURCE|HGNC_ID|LoF|LoF_filter|LoF_flags|LoF_info'.split("|")
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

def query_tabix(ch, s, e, opt=""):
    os.system(f"tabix -h {CADD_FILE} {ch}:{s}-{e} > /tmp/{ch}:{s}-{e}.tsv")
    return pd.read_table(f"/tmp/{ch}:{s}-{e}.tsv", header=1)

def get_cadd(ch,pos,ref,alt):
    # Only supports SNVs at this time
    if not (len(ref) == 1 and len(alt)== 1): 
        return [""] * len(cadd_anno)

    # Query tabix
    cadd_table = query_tabix(ch, pos, int(pos))
    cadd_table = cadd_table[cadd_table["Alt"] == alt][cadd_anno]
    if not len(cadd_table): 
        return [""] * len(cadd_anno)
    
    # TODO pick the most impactful (min or max) for each variant
    cadd = cadd_table.iloc[0].array.astype(str).tolist()
    
    # Replace "nan" with ""
    svi = np.where(cadd_table.columns == "SIFTval")[0][0]
    if cadd[svi] == "nan": cadd[svi] = ""
    pvi = np.where(cadd_table.columns == "PolyPhenVal")[0][0]
    if cadd[pvi] == "nan": cadd[pvi] = "" 

    # Replace SIFT with categoriacal
    sci = np.where(cadd_table.columns == "SIFTcat")[0][0]
    sc = to_cat(cadd[sci], "SIFTcat").values.reshape(-1).astype(str).tolist()
    cadd[sci] = sc[0]
    for e in range(len(sc[1:])): cadd.insert(sci+e,sc[e])

    # Replace PolyPhen with categorical    
    pci = np.where(cadd_table.columns == "PolyPhenCat")[0][0] + len(sc[1:])
    pc = to_cat(cadd[pci], "PolyPhenCat").values.reshape(-1).astype(str).tolist()
    cadd[pci]=pc[0]
    for e in range(len(pc[1:])): cadd.insert(pci+e,pc[e])
    return cadd

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
    
def process_line(l):
    l = l.strip("\n").split()
    rare_idx = np.where(np.isin(l, ["0/1", "1/0", "1/1", "0|1", "1|0", "1|1"]))[0]
    rare_ids = [cols[i] for i in rare_idx]
    af = l[7].split("AFR_AF=")[1].split(";")[0]

    vep_cat_df = get_vep(l)
    if vep_cat_df is None: return None

    ch,pos,ref,alt = l[chidx], l[pidx], l[refidx], l[altidx]
    cadd = get_cadd(ch,pos,ref,alt)
    
    out_lines = []

    # Create annotation line for each id, gene pair
    for rid in set(rare_ids).intersection(inds_keep):
        for gene in vep_cat_df.index:
            
            # Filter out indices with too many outliers
            if not rid in inds_keep: print(rid)
            # TODO For now I'm skipping genes that aren't in expression data
            # but I think its okay to have some blank values
            if (gene not in eout_df): continue

            out_line = [rid, gene, ch, pos, af]
            out_line.extend(vep_cat_df.loc[gene].array.astype(str).tolist())

            out_line.extend(cadd)
            eout = eout_df.loc[rid, gene]
            out_line.append(str(eout))
            out_lines.append("\t".join(out_line) + "\n")

    genes = vep_cat_df.index.array
    return "".join(out_lines), genes, rare_ids


def get_header():
    header = ["SubjectID", "GeneName", "Chromosome", "Position", "AF"]
    header.extend(cats["Consequence"])
    header.extend(cats["LoF"])

    cadd_anno_cols = []
    for c in cadd_anno:
        if c[-3:].lower() ==  "cat":
            cadd_anno_cols.extend([f"{c[:-3]}_{cc}" for cc in cats[c]])
        else: cadd_anno_cols.append(c)
    header.extend(cadd_anno_cols)
    header.append("eOutlier") 
    return "\t".join(header) + "\n"

vcf = open(vcf_file, 'r')
cols = vcf.readline()[1:-1].split()

chidx = np.where(np.array(cols) == "CHROM")[0][0]
pidx = np.where(np.array(cols) == "POS")[0][0]
refidx = np.where(np.array(cols) == "REF")[0][0]
altidx = np.where(np.array(cols) == "ALT")[0][0]

lidx = np.where(np.array(vep_fields) == "LoF")[0][0]
vidx = np.where(np.array(vep_fields) == "Consequence")[0][0]
gidx = np.where(np.array(vep_fields) == "Gene")[0][0]

out = open(f"{data_dir}/{tsv_out}", 'w')
gout = open(f"{data_dir}/{gene_outliers}", 'w')
print(f"{data_dir}/{tsv_out}")
print(f"{data_dir}/{gene_outliers}")

# Expression outlier scores from residuals file
eout_file = f"/oak/stanford/groups/smontgom/erobb/data/watershed/{pop}_exprResiduals.tsv"
eout_keep_file = f"{data_dir}/ids_outlier_filtered_{pop}_t3f75.txt"
eout_df = pd.read_table(eout_file, sep="\t", index_col=0).T
# drop transcript information
gene_names = [c.split(".")[0] for c in eout_df.columns]
assert(len(np.unique(gene_names)) == len(eout_df.columns))
eout_df.columns = gene_names
inds_keep = pd.read_table(eout_keep_file, sep=" ", index_col=0, header=None).T
#eout_df = eout_df.loc[eout_df.index.intersection(inds_keep.columns)]

    
header = get_header()
out.write(header)

def do_work(in_queue, out_list):
    t1  = time.time()
    while True:
        if not in_queue.empty():
            item = in_queue.get(timeout=1)
            line_no, line = item
            # exit signal 
            if line == None:
                return
            result = process_line(line)
            out_list.append((line_no, result))
        # timeout
        elif time.time() - t1 > 1: return

if __name__ == "__main__":

    os.system(f"wc -l {vcf_file} > /tmp/wc_{pop}")
    n = int(open(f"/tmp/wc_{pop}", "r").readlines()[0].split()[0])
    num_workers = args.nw
    lines_ps = args.buf

    manager = Manager()
    results = manager.list()
    work = manager.Queue(num_workers)

    pool = []
    for i in range(num_workers):
        p = Process(target=do_work, args=(work, results))
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
                p = Process(target=do_work, args=(work, results))
                p.start()
                pool.append(p)

    out.close()  
    gout.close()



