#!/bin/bash


data_dir="${oak}/data/watershed"
gencodeprefix="gencode.v43.chr_patch_hapl_scaff.annotation"
gencode="gencode.v43.chr_patch_hapl_scaff.annotation.gtf.gz"

# for gencode annotation get lincRNA and prootein-coding exons, but also pad the internal exons by 5 base pairs
python pad.gtf.exons.py $gencode | sort -k1,1 -k2,2n | uniq > ${gencodeprefix}_coding.lincRNA_padded.bed


