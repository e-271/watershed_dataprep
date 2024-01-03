#!/usr/bin/env Rscript
library(dplyr)

args = commandArgs(trailingOnly=TRUE)
tsv=args[1]
rename=read.csv(args[2],header=FALSE,row.names=1,na="")

df=read.table(tsv, fill=TRUE, header=1, sep="\t")

drop = row.names(rename %>% filter(is.na(V2)))
rename = rename %>% filter(!is.na(V2))

# rename
df <- df %>% rename_with(~rename$V2, all_of(row.names(rename)))
# drop
df <- df %>% select(-all_of(drop))

test <- df %>% filter(is.na(N2pair))
train <- df %>% filter(!is.na(N2pair))

pfx = tools::file_path_sans_ext(tsv)
write.table(df,sprintf("%s.format.tsv", pfx),
            quote=FALSE,row.names=FALSE, na="NA", sep="\t")
write.table(train,sprintf("%s.format.train.tsv", pfx),
            quote=FALSE,row.names=FALSE, na="NA", sep="\t")
write.table(test,sprintf("%s.format.test.tsv", pfx),
            quote=FALSE,row.names=FALSE, na="NA", sep="\t")

