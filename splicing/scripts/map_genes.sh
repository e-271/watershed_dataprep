
cts=$1
gc=$2
filt_juncs=$3
intersect_out=$4

# Convert Leafcutter output to .junc format for bedtools
tail -n +2 $cts | cut -d' ' -f 1 | \
    sed 's/:/\t/g' | awk -F'\t' '{print $1 FS $2 FS $3-1 FS $4}'  > $filt_juncs

# Intersect Leafcutter juntions with genocde gene annotations
bedtools intersect -f 1 -a $filt_juncs -b $gc -loj | cut -f 1-4,9 > $intersect_out

# Filter out junctions with no gene annotation
cat $intersect_out | cut -f 4,5 | sort | uniq | awk '{if ($2 != ".") print $0}' 

