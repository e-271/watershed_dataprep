#!/bin/bash

cadd=$1
cadd_indel=$2
cadd_cols=$3
vcf=$4


# Tell bcftools what columns we want from CADD file, based on input.cols (assumes 1 column per line)
tabix -h $cadd 1:1-1 | tail -n 1 | sed 's/#//' | sed 's/\s/,/g' > all_columns.tmp
cols=$(awk 'NR==FNR{a[$0];next}{for(i=1; i<=NF; i++) {if($i in a){printf $i","} else if($i == "Chrom" || $i == "Pos" || $i == "Ref" || $i == "Alt"){printf $i","}else{printf "-,"}}}' $cadd_cols FS="," all_columns.tmp)

# Create VCF header info lines for relevant CADD annotations
merge_logic=""
while read anno; do
    echo "##INFO=<ID=${anno},Number=.,Type=String,Description=\"CADD annotation\">" >> header.tmp
    merge_logic="${merge_logic}${anno}:max,"
done < $cadd_cols

for l in {1..22} X Y; do echo chr$l $l; done > rename_chr.tmp
bcftools annotate --rename-chrs rename_chr.tmp $vcf \
| bcftools annotate -a $cadd -c ${cols} -h header.tmp --merge-logic ${merge_logic} --min-overlap 1 \
| bcftools annotate -a $cadd_indel -c ${cols} -h header.tmp --merge-logic ${merge_logic} --min-overlap 1 

rm rename_chr.tmp header.tmp all_columns.tmp


