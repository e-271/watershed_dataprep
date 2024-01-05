#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
tsv=args[1]
config=args[2]

drop_column <- function(df, col_name) {
    j=which(colnames(df) == col_name)
    return(cbind(df[1:(j-1)], df[(j+1):dim(df)[2]]))
}


# Convert string column to binary categorical vector.
# Multiple categories may be represented in the same entry, eg: cat1,cat2,cat3
to_cat = function(df, col_name, drop_cols) {
    df[,col] = as.character(df[,col])
    col_split=strsplit(df[,col_name], ',')
    cats=unique(unlist(col_split))
    cats=cats[!is.na(cats) & !("" == cats)]
    if (length(drop_cols)) { cats=cats[!(cats %in% drop_cols)] }
    if (length(cats) < 1) {
        warning(sprintf("No unique valid value in column %s. Dropping.", col_name))
        df_cat <- drop_column(df, col_name)
        if (col_name == "SIFTcat") { df_cat <- drop_column(df_cat, "SIFTval") } 
        if (col_name == "PolyPhenCat") { df_cat <- drop_column(df_cat, "PolyPhenVal") } 
        return(df_cat)
    }
    j=which(colnames(df) == col_name)
    col_binary=sapply(cats, function(x) sapply(col_split, function(y) x %in% y))
    colnames(col_binary)=sapply(cats, function(x) paste0(col_name,"_",x))
    df_cat=cbind(df[1:(j-1)], col_binary*1, df[(j+1):dim(df)[2]])
    return(df_cat)
}

df=read.table(tsv, header=1,sep="\t",na=c("", "NA"))
categ=read.table(config,row.names=1,header=1,fill=TRUE)
for (col in row.names(categ)) {
  drop_cols = unlist(strsplit(categ[col,],","))
  df=to_cat(df,col,drop_cols)
}
write.table(df,"",quote=FALSE,row.names=FALSE,na="",sep="\t")

