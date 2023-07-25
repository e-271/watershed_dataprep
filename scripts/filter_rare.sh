#!/bin/bash

DATA_DIR='/oak/stanford/groups/smontgom/erobb/data/vep'
POP=$1
OUT="AF.all.${POP}.hg38a.ID.ba.VEP.vcf"
OUT_RARE="AF.all.${POP}.hg38a.ID.ba.VEP.rare.vcf"
head -n 100 ${DATA_DIR}/${OUT}  | awk '$0 ~ /^#[a-z,A-Z,0-9]/' > ${DATA_DIR}/${OUT_RARE} 
cat ${DATA_DIR}/${OUT} | \
awk '{if ($0 ~ /^[^#]/) {print $0}}' | \
awk '{n = split($8, a, ";"); af = substr(a[7], 8); if (af <= 0.01 && af > 0) {print $0} }'  >> ${DATA_DIR}/${OUT_RARE}
#awk '{n = split($8, a, \"=\"); if (a[2] < 0.01 && a[2] > 0) {print $0}}' >> ${DATA_DIR}/${OUT_RARE}"

