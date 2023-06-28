#!/bin/bash

data_dir="/oak/stanford/groups/smontgom/erobb/data"

cat "${data_dir}/gencode/gencode.v43.chr_patch_hapl_scaff.annotation.gtf" |  awk -F"\t" '{if ($3 == "transcript") {print $0} }' | awk -F"\t" '{n=split($9,a,";"); b=substr(a[3],12); if (b=="\"protein_coding\"" || b=="\"lincRNA\"" ) {print $0} }' > "${data_dir}/gencode/gencode.v43.chr_patch_hapl_scaff.annotation.transcripts.protein_lincRNA.gtf"

cat "${data_dir}/gencode/gencode.v43.chr_patch_hapl_scaff.annotation.gtf" |  awk -F"\t" '{if ($3 == "transcript") {print $0} }' | awk -F"\t" '{n=split($9,a,";"); b=substr(a[3],12); if (b=="\"protein_coding\"" || b=="\"lincRNA\"" ) {print substr(a[1],10)} }' |  sed 's/"//g' | sort | uniq >  "${data_dir}/gencode/gencode.v43.gene_ids.protein_lincRNA.txt"

#cat "${data_dir}/gencode/gencode.v43.chr_patch_hapl_scaff.annotation.gtf" |  awk -F"\t" '{if ($3 == "transcript") {print $0} }' > "${data_dir}/gencode/gencode.v43.chr_patch_hapl_scaff.annotation.transcripts.gtf"


