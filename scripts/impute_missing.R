#!/usr/bin/env Rscript
library(dplyr)

args = commandArgs(trailingOnly=TRUE)
tsv=args[1]
impute=read.csv(args[2],header=FALSE, row.names=1)
nout=as.numeric(args[3])

df=read.table(tsv, header=1, sep="\t",na=c("", "NA"))

# Replace missing values with imputation values
impute_cols=intersect(colnames(df), row.names(impute))

missing=setdiff(row.names(impute),impute_cols)
if (length(missing)) {
warning(paste(c("Columns from config/impute not found in data:", missing), collapse=" "))
}
missing=setdiff(colnames(df)[3:(length(colnames(df))-2-nout)],impute_cols)
if (length(missing)) {
warning(paste(c("Columns from data not found in config/impute:", missing), collapse=" "))
}

for (col in impute_cols){
  df[is.na(df[,col]),col] = impute[col,]
}


write.table(df,"",quote=FALSE,row.names=FALSE, na="NA", sep="\t")
