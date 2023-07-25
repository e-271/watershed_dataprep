#!/bin/bash

DATA_DIR=/oak/stanford/groups/smontgom/erobb/data

cat ${DATA_DIR}/gencode/gencode.v43.chr_patch_hapl_scaff.annotation.gtf \
| awk '{if($3 == "gene") {print $0}}' \
> ${DATA_DIR}/gencode/gencode.v43.chr_patch_hapl_scaff.annotation.genes.tsv


cat ${DATA_DIR}/gencode/gencode.v43.chr_patch_hapl_scaff.annotation.gtf \
| awk '{if($3 == "gene") {print $0}}' \
| cut -f9,1,4,5 \
| awk -F"\t" '{n = split($4, a, ";"); print substr(a[1],9) FS $1 FS $2 FS $3}' \
| sed 's/"//g' \
> ${DATA_DIR}/gencode/gencode.v43.gene_pos.tsv


