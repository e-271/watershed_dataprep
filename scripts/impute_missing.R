#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
tsv=args[1]
impute=read.csv(args[2],header=FALSE, row.names=1)

df=read.table(tsv, fill=TRUE, header=1, sep="\t")

# Replace missing values with imputation values
for (col in row.names(impute)[1:2]) {
  if (!(col %in% colnames(df))) {
        warning(paste(col,"not found."))
        next
  }
  df[is.na(df[,col]),col] = impute[col,]
}

write.table(df,"",quote=FALSE,row.names=FALSE, na="", sep="\t")

