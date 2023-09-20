#!/bin/bash

vcf=$1
fields=$2

format=$(cat $fields | awk -F, '{ORS=""; print "\\t" "%" $1}')
format="[%SAMPLE\\t%INFO/Gene${format}\n]"

header=$(cat $fields | awk -F, 'BEGIN{printf "Sample\tGene"}{ printf "\t" $1}')
echo "$header"

inc="GT=\"alt\""
bcftools query -i "$inc" -r "1"  -f "$format" $vcf


