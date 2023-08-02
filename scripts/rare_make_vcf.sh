#!/bin/bash

POP=$1
DATA_DIR='/oak/stanford/groups/smontgom/erobb/data/vcf/30x'
AF="all.${POP}.AFGR.30x.frq"
VCF="all.${POP}.30x.ID.vcf"
OUT_RARE="all.${POP}.30x.ID.rare_afgr.vcf"


rm  ${DATA_DIR}/frq/${AF}.vcf

head -n 5000 ${DATA_DIR}/${VCF}  | awk '$0 ~ /^##*[a-z,A-Z,0-9]/' > ${DATA_DIR}/frq/${AF}.vcf

cat "${DATA_DIR}/frq/${AF}" | \
awk '{if ($3 <= 0.01 && $2 > 0) {split($1, a, "_"); print a[1] "\t" a[2] "\t" a[3] "\t" a[4] "\t\t\tAFGR_AF="$3} }' \
>> ${DATA_DIR}/frq/${AF}.vcf 

#awk '{if ($2 <= 0.01 && $2 > 0) {split($1, a, "_"); print a[1] "\t" a[2] "\t" a[3] "\t" a[4] "\t\t\tAF="$2} }' \

head -n 5000 ${DATA_DIR}/${VCF}  | awk '$0 ~ /^#[a-z,A-Z,0-9]/' > ${DATA_DIR}/${OUT_RARE}
bedtools intersect -a ${DATA_DIR}/${VCF} -b ${DATA_DIR}/frq/${AF}.vcf -f 1 -F 1 -wa >> ${DATA_DIR}/${OUT_RARE}

