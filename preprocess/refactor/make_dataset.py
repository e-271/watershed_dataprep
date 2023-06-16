import numpy as np
import os
import pandas as pd
import numpy as np
import argparse

import numpy as np
import os
import pandas as pd
import numpy as np


argParser = argparse.ArgumentParser()
argParser.add_argument("--pop", default="AFR", type=str)
argParser.add_argument("--vep_dir", default='/oak/stanford/groups/smontgom/erobb/data/vep', type=str)
argParser.add_argument("--data_dir", default='/oak/stanford/groups/smontgom/erobb/data/watershed', type=str)
argParser.add_argument("--cadd", default="/oak/stanford/groups/smontgom/erobb/data/watershed/whole_genome_SNVs_inclAnno.tsv.gz", type=str)

args = argParser.parse_args()
pop = args.pop
vep_dir, data_dir = args.vep_dir, args.data_dir
CADD_FILE = args.cadd


vcf_in = f"AF.all.{pop}.hg38a.ID.ba.VEP.rare.vcf"
tsv_out = f'AF.all.{pop}.hg38a.ID.ba.VEP.rare.ws.tsv'
gene_outliers = f'gene_outliers_{pop}.tsv'
vcf_file = f"{vep_dir}/{vcf_in}"
eout_file = f"{data_dir}/eOutlier_scores_{pop}_t3.txt"
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
    if var == 'nan' and cat == 'PolyPhenCat': var = 'unknown'
    return np.isin(cats[cat],var).astype(int).astype(str).tolist()

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
    sc = to_cat(cadd[sci], "SIFTcat")
    cadd[sci] = sc[0]
    for e in range(len(sc[1:])): cadd.insert(sci+e,sc[e])

    # Replace PolyPhen with categorical    
    pci = np.where(cadd_table.columns == "PolyPhenCat")[0][0] + len(sc[1:])
    pc = to_cat(cadd[pci], "PolyPhenCat")
    cadd[pci]=pc[0]
    for e in range(len(pc[1:])): cadd.insert(pci+e,pc[e])
    return cadd

def process_line(l):
    l = l.strip("\n").split()
    rare_idx = np.where(np.isin(l, ["0/1", "1/0", "1/1", "0|1", "1|0", "1|1"]))[0]
    rare_ids = [cols[i] for i in rare_idx]

    if 'CSQ=' not in l[7]: return None
    af = l[7].split("AFR_AF=")[1].split(";")[0]
    vep_infos = l[7].split('CSQ=')[1].split(',')
    vep_infos = [v.split('|') for v in vep_infos]
    for v in vep_infos:
        # Get VEP variant and LoF
        var = to_cat(v[vidx].split("&"), "Consequence")
        lof = to_cat(v[lidx], "LoF")
        gene = v[gidx]

    ch,pos,ref,alt = l[chidx], l[pidx], l[refidx], l[altidx]
    cadd = get_cadd(ch,pos,ref,alt)
    
    out_lines = []
    # Skip individuals with > threshold outlier genes
    for rid in set(rare_ids).intersection(inds_keep):
        if not rid in inds_keep: print(rid)
        # TODO For now I'm skipping genes that aren't in expression data
        # but I think its okay to have some blank values
        if (gene not in eout_df): continue
        out_line = [rid, gene, ch, pos, af]
        out_line.extend(var)
        out_line.extend(lof)
        out_line.extend(cadd)
        eout = eout_df.loc[rid, gene]
        out_line.append(str(eout))
        out_lines.append("\t".join(out_line) + "\n")
    
    return "".join(out_lines), gene, rare_ids


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

eout_df = pd.read_table(eout_file, sep=" ")
# drop transcript information
gene_names = [c.split(".")[0] for c in eout_df.columns]
assert(len(np.unique(gene_names)) == len(eout_df.columns))
eout_df.columns = gene_names
inds_keep = pd.read_table(eout_keep_file, sep=" ", index_col=0, header=None).T
    
header = get_header()
out.write(header)

for l in vcf:
    out_line, gene, rid = process_line(l)
    if out_line is not None and out_line != "":
        #print(out_line)
        out.write(out_line)
        gout.write("\t".join([gene] + rid) + "\n")
        
out.close()  
gout.close()



