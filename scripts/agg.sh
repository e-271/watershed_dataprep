#!/bin/bash

vcf_dir=$1
fields=$2
gencode=$3
proc=$4

format=$(cat $fields | awk -F, '{ORS=""; print "%" $1 "\t"}')
format="[${format}\n]"

header=$(cat $fields | awk -F, 'BEGIN{printf "Sample\tGene"}{ printf "\t" $1}')
echo "$header"

query () {
    vcf=$1
    gene=$2
    chr=$3
    start=$4
    end=$5

    source scripts/agg_functions.sh
    inc="INFO/Gene=\"${gene}\"&&GT=\"alt\""
    id=$(basename $vcf .vcf.gz)

    q=$(bcftools query -i "$inc" -r $chr:$((start-10000))-$((end+10000)) -f "$format" $vcf)
    if [ "${#q}" == 0 ]; then return; fi

    i=1
    line="${id}\t${gene}"
    while IFS=',' read field func; do
        agg=$(echo "$q" | cut -f$i | split | clean_missing | eval $func)
        line="${line}\t${agg}"
        ((i++))
    done < $fields
    printf "${line}\n"

}
export -f query
export format fields

# TODO it may be much faster not to split by ID, to reduce the number of queries
while IFS=' ' read -r gene chr start end; do 
    gene=$(echo $gene | sed 's/\.[0-9]*//') # remove .[0-9]* from gene ID to match VEP
    # Sequential
    #for vcf in $vcf_dir/*.vcf.gz; do
    #    query $vcf $gene $chr $start $end
    #done
    # Parallel
    ls $vcf_dir/*.vcf.gz | xargs -P $proc -n 1 sh -c "query \"\$1\" $gene $chr $start $end" {}
done < $gencode

