#!/bin/bash

# Download human_ancestor.fa.gz, phylocsf_gerp.sql from 
# https://github.com/konradjk/loftee


VEP_DIR="/oak/stanford/groups/smontgom/erobb/data/vep"
VCF_DIR="/oak/stanford/groups/smontgom/pgoddard/africa/data/genotypes/vcf"
POP=$1
echo $1

cd ~/ensembl/ensembl-vep
./vep \
-i $VCF_DIR/AF.all.${POP}.hg38a.ID.ba.vcf.gz \
-o $VEP_DIR/AF.all.${POP}.hg38a.ID.ba.VEP.vcf \
--verbose \
--vcf \
--fork 64 \
--no_stats \
--force_overwrite \
--offline \
--dir_cache $VEP_DIR \
--plugin LoF,\
human_ancestor_fa:$VEP_DIR/human_ancestor.fa.gz,\
loftee_path:/home/erobb/.vep/Plugins,\
conservation_file:$VEP_DIR/phylocsf_gerp.sql



