#!/bin/bash

TXID=9
PVAL=14
PHET=8

data_dir=/oak/stanford/groups/smontgom/erobb/data

POP=$1


# Copy header
rm "${data_dir}/eqtls/sorted.byGene.dist.hwe.af.${POP}.eQTL.nominal.hg38a.txt"
head -n 1  "${data_dir}/eqtls/sorted.dist.hwe.af.${POP}.eQTL.nominal.hg38a.txt" \
>  "${data_dir}/eqtls/sorted.byGene.dist.hwe.af.${POP}.eQTL.nominal.hg38a.txt"

echo "Sorting..."
# Sort by gene & pvalue
cat  "${data_dir}/eqtls/sorted.dist.hwe.af.${POP}.eQTL.nominal.hg38a.txt" | \
tail -n +2 | \
sort -S 50% --parallel=32 -b -k${TXID},${TXID} -k${PHET}g \
>>  "${data_dir}/eqtls/sorted.byGene.dist.hwe.af.${POP}.eQTL.nominal.hg38a.txt"

# Copy header
rm  "${data_dir}/eqtls/sorted.byGene.dist.hwe.af.${POP}.top1.eQTL.nominal.hg38a.txt"
head -n 1  "${data_dir}/eqtls/sorted.dist.hwe.af.${POP}.eQTL.nominal.hg38a.txt" \
>  "${data_dir}/eqtls/sorted.byGene.dist.hwe.af.${POP}.top1.eQTL.nominal.hg38a.txt"

echo "Extracting top eQTL..."
# Take lowest pvalue for each gene
cat  "${data_dir}/eqtls/sorted.byGene.dist.hwe.af.${POP}.eQTL.nominal.hg38a.txt" | \
tail -n +2 | \
awk '{if (C != $9) {C=$9; print $0}}' \
>>  "${data_dir}/eqtls/sorted.byGene.dist.hwe.af.${POP}.top1.eQTL.nominal.hg38a.txt"

echo "Extracting ID/variant for top eQTLs..."
# Copy header
rm  "${data_dir}/eqtls/AF.all.${POP}.hg38aID.eQTLs.ba.vcf"
head -n 100 "${data_dir}/vcf/AF.all.${POP}.hg38a.ID.ba.vcf" | awk '{if ($1 == "#CHROM"){print "GENE" FS $0}}' | sed s/#//  >  "${data_dir}/eqtls/AF.all.${POP}.hg38aID.eQTLs.ba.vcf"
# Query vcf for each top1 eQTL
cat  "${data_dir}/eqtls/sorted.byGene.dist.hwe.af.${POP}.top1.eQTL.nominal.hg38a.txt" | awk -v pop="${POP}" -v data="${data_dir}" '{"tabix " data "/vcf/AF.all." pop ".hg38a.ID.ba.vcf.gz " $1 ":" $2 "-" $3 | getline n; if(n){print $9 FS n}}' >> "${data_dir}/eqtls/AF.all.${POP}.hg38aID.eQTLs.ba.vcf"




