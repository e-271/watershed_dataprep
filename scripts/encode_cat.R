#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
tsv=args[1]
config=args[2]

# Convert string column to binary categorical vector.
# Multiple categories may be represented in the same entry, eg: cat1,cat2,cat3
to_cat = function(df, col_name, drop_cols) {
    cons=strsplit(df[,col_name], ',')
    f=unique(unlist(cons))
    f=f[!is.na(f) & !("" == f)]
    if (length(drop_cols)) {
        f=f[!(f %in% drop_cols)]
    }
    bin=sapply(f, function(x) sapply(cons, function(y) x %in% y))
    colnames(bin)=sapply(f, function(x) paste0(col_name,"_",x))
    j=which(colnames(df) == col_name)
    df2=cbind(df[1:(j-1)], bin*1, df[(j+1):dim(df)[2]])
    return(df2)
}

df=read.table(tsv, fill=TRUE, header=1,sep="\t",na=c("", "NA"))
categ=read.table(config,row.names=1,header=1,fill=TRUE)
for (col in row.names(categ)) {
  drop_cols = unlist(strsplit(categ[col,],","))
  df=to_cat(df,col,drop_cols)
}
write.table(df,"",quote=FALSE,row.names=FALSE,na="",sep="\t")

