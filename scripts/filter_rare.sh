#!/bin/bash

POP=$1
DATA_DIR='/oak/stanford/groups/smontgom/erobb/data/vcf'
#OUT="AF.all.${POP}.hg38a.ID.ba.vcf"
#OUT_RARE="AF.all.${POP}.hg38a.ID.ba.rare.vcf"
OUT="AF.all.${POP}.hg38a.ID.matt.vcf"
OUT_RARE="AF.all.${POP}.hg38a.ID.matt.rare.vcf"


#rm ${DATA_DIR}/${OUT_RARE}
head -n 100 ${DATA_DIR}/${OUT}  | awk '$0 ~ /^#[a-z,A-Z,0-9]/' > ${DATA_DIR}/${OUT_RARE} 
cat ${DATA_DIR}/${OUT} | \
awk '{if ($0 ~ /^[^#]/) {print $0}}' | \
awk '{n = split($8, a, ";"); af = substr(a[4], 4); if (af <= 0.01 && af > 0) {print $0} }' >> ${DATA_DIR}/${OUT_RARE}

echo ${DATA_DIR}/${OUT_RARE}

#awk '{n = split($8, a, \"=\"); if (a[2] < 0.01 && a[2] > 0) {print $0}}' >> ${DATA_DIR}/${OUT_RARE}"

