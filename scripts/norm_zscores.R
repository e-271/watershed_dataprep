#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
tsv=args[1]
zscores=unlist(strsplit(args[2],','))

df=read.table(tsv, fill=TRUE, header=1)

# Replace z-scores with norm (so they are comparable to p-values)
for (zs in zscores) {
  df[,zs] = dnorm(df[,zs])
}

write.table(df,"",quote=FALSE,row.names=FALSE)

