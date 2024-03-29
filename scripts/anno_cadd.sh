#!/bin/bash

cadd=$1
cadd_cols=$2
vcf=$3

# Tell bcftools what columns we want from CADD file, based on input.cols (assumes 1 column per line)
tabix -h $cadd 1:1-1 | tail -n 1 | sed 's/#//' | sed 's/\s/,/g' > ${vcf}.all_columns.tmp
# Mark CADD columns to keep and discard
cols=$(
awk 'NR==FNR{a[$1];next}
{for(i=1; i<=NF; i++) {
    if($i in a || ($i == "Chrom" || $i == "Pos" || $i == "Ref" || $i == "Alt")){
        printf $i","
    }
    else{printf "-,"}
}
}' FS=',' $cadd_cols FS="," ${vcf}.all_columns.tmp)

# Create VCF header info lines for relevant CADD annotations
merge_logic=""
while IFS=',' read -r anno merge; do
    echo "##INFO=<ID=${anno},Number=.,Type=String,Description=\"CADD annotation\">" >>  ${vcf}.header.tmp
    merge_logic="${merge_logic}${anno}:${merge},"
done < $cadd_cols

bcftools annotate $vcf -a $cadd -c ${cols} -h ${vcf}.header.tmp --merge-logic ${merge_logic} --min-overlap 1 

rm  ${vcf}.header.tmp  ${vcf}.all_columns.tmp

