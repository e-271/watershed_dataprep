library(dplyr)

args = commandArgs(trailingOnly=TRUE)
counts=args[1]
gencode=args[2]
min_ct=as.numeric(args[3])

cts <- read.csv(counts, header=T, sep=",")
gc <- read.table(gencode, col.names=c("Chrom", "Start", "End", "Strand", "Gene", "Name", "Type"))

# filter to protein-coding / lincRNA genes
d1=dim(cts)[1]
gc <- gc %>% filter(Type == "lincRNA" | Type == "protein_coding")
cts <- cts %>% filter(Gene %in% gc$Gene)
d2=dim(cts)[1]
warning(sprintf("filtered out %d/%d non protein-coding/lincRNA genes", d1-d2, d1))

# drop genes with few reads
d1=dim(cts)[1]
cts <- cts %>% filter(rowSums(across(where(is.numeric))) > min_ct )
d2=dim(cts)[1]
warning(sprintf("filtered out %d/%d genes with low read counts", d1-d2, d1))

write.csv(cts,"",row.names=F)

