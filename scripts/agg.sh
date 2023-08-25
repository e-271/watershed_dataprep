#!/bin/bash

vcf_dir=$1
fields=$2
gencode=$3

alias max="sort -n | tail -n 1"
alias min="sort -n | head -n 1"
alias unique="sort -u | tr '\n' ',' | sed 's/,$//'"
alias append="tr '\n' ',' | sed 's/,$//'"

alias split="sed 's/[,&]/\n/g'"
alias clean_missing="sed 's/^\.$//g' | sed 's/^NA$//g'"

format=$(cat $fields | awk -F, '{ORS=""; print "%" $1 "\t"}')
format="[${format}\n]"

header=$(cat $fields | awk -F, 'BEGIN{printf "Sample\tGene"}{ printf "\t" $1}')
echo "$header"

while IFS=' ' read -r gene chr start end; do 
    gene=$(echo $gene | sed 's/\.[0-9]*//') # remove .[0-9]* from gene ID to match VEP
    chr=$(echo $chr | sed 's/chr//') # remove chr to match our vcf (done in an earlier snakemake step)
    for vcf in $vcf_dir/*.vcf.gz; do
        inc="INFO/Gene=\"${gene}\"&&GT=\"alt\""
        id=$(basename $vcf .vcf.gz)
        
        q=$(bcftools query -i "$inc" -r $chr:$((start-10000))-$((end+10000)) -f "$format" $vcf)
        if [ "${#q}" == 0 ]; then continue; fi

        i=1
        line="${id}\t${gene}"
        while IFS=',' read field func; do
            agg=$(echo "$q" | cut -f$i | split | clean_missing | eval $func)
            line="${line}\t${agg}"
            ((i++))
        done < $fields
        printf "${line}\n"

        # This is slower unfortunately, doing 20 bcftools query vs. 20 cut -f is probably inefficient.
        #while IFS=',' read field func; do
        #    agg=$(bcftools query -i "$inc" -r $chr:$((start-10000))-$((end+10000)) -f "%$field\\n" $vcf | eval $func )
        #    line="${line}\t${agg}"
        #done < $fields
        #printf "${line}\n"
    done
done < $gencode

