input=$1
min_reads=$2

# The cluster boundaries list [start, end] but junc files list [start, end-1] so we subtract 1 to match.
tail -n +2 $input | \
            awk -F' ' -v min=$2 '{s=0; {for(i=2;i<=NF;i++){ if ($i > s) {s = $i} }}; if(s>min){print $1}} ' | \
            sed 's/:/\t/g' | \
            awk -F'\t' '{print $1 FS $2 FS $3-1}'

