#!/bin/bash

input=$1
name=$2

echo -e "Sample\tGene\t${name}"

cat $input | \
awk -F'\t' 'BEGIN{OFS="\t"}{
sub(/\.[0-9][0-9]*/, "", $1); 
print $0}' | \
awk  -F'\t' 'BEGIN{OFS="\t"}{
if(NR==1)
    {split($0,a,FS)}
else{
    for(i=2; i<NF; i++){
        print a[i-1] OFS $1 OFS $i }
    }
}'


