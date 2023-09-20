library(tidyr)
library(dplyr)

args = commandArgs(trailingOnly=TRUE)
 
tsv=args[1]
keepgenes = unlist(read.table(args[2], header=F))
zthresh=args[3]
nout_thresh=args[4]

df=read.table(tsv, header=T, sep="\t", na="")

# Genes with at least 1 outlier
has_outliers <- df %>% group_by(Gene) %>% 
    filter(abs(eOutliers) > zthresh) %>% 
    summarise()
has_outliers <- unlist(has_outliers)

# Samples with <threshold eOutliers
under_outlier_thresh <- df %>% group_by(Sample) %>%
    summarise(noutliers=sum(abs(eOutliers) > zthresh)) %>%
    filter(noutliers < nout_thresh)
under_outlier_thresh = under_outlier_thresh$Sample

df_keep <- df %>% filter(Gene %in% keepgenes & Gene %in% has_outliers) %>% filter(Sample %in% under_outlier_thresh)

write.table(df_keep,"",quote=FALSE,row.names=FALSE,na="",sep="\t")

