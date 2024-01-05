#!/usr/bin/env Rscript
library(dplyr)

args = commandArgs(trailingOnly=TRUE)
tsv=args[1]
df=read.table(tsv, fill=TRUE, header=1, sep="\t")

test <- df %>% filter(!is.na(N2pair))
train <- df %>% filter(is.na(N2pair))

pfx = tools::file_path_sans_ext(tsv)
write.table(train,sprintf("%s.train.tsv", pfx),
            quote=FALSE,row.names=FALSE, na="NA", sep="\t")
write.table(test,sprintf("%s.test.tsv", pfx),
            quote=FALSE,row.names=FALSE, na="NA", sep="\t")

