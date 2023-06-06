#!/bin/bash

cat /oak/stanford/groups/smontgom/erobb/data/watershed/gencode.v43.chr_patch_hapl_scaff.annotation.gtf |  awk -F"\t" '{if ($3 == "exon") {print $0} }' | awk -F"\t" '{n=split($9,a,";"); b=substr(a[3],12); if (b=="\"protein_coding\"" || b=="\"lincRNA\"" ) {print $0} }' > /oak/stanford/groups/smontgom/erobb/data/watershed/gencode.v43.chr_patch_hapl_scaff.annotation.exons.protein_lincRNA.gtf

#TODO pad 
 cat /oak/stanford/groups/smontgom/erobb/data/watershed/gencode.v43.chr_patch_hapl_scaff.annotation.exons.protein_lincRNA.gtf | awk -F"\t" '{n=split($9,a,";"); s=""; for (i in a) {split(a[i],b," "); s = s FS b[2]}; print s }' > /oak/stanford/groups/smontgom/erobb/data/watershed/gencode.v43.chr_patch_hapl_scaff.annotation.exons.protein_lincRNA.bed


