#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
tsv=args[1]
cols=unlist(strsplit(args[2],','))
tsv_out=args[3]

#tsv="data/watershed/30x.LWK.filt.ref_af.rare.CADD.VEP_split.id_split.agg.eOutliers.pairlabel.tsv"
#cols=unlist(strsplit("Consequence,LoF,SIFTcat,PolyPhenCat",','))
#tsv_out="data/watershed/30x.LWK.filt.ref_af.rare.CADD.VEP_split.id_split.agg.eOutliers.pairlabel.cat.tsv"

# Convert string column to binary categorical vector.
# Multiple categories may be represented in the same entry, eg: cat1,cat2,cat3
to_cat = function(df, col_name, drop_cols) {
    cons=strsplit(df[,col_name], ',')
    f=unique(unlist(cons))
    f[!(f %in% drop_cols) & !is.na(f)]
    cats = !is.na(f) & !(f==".")
    bin=sapply(f, function(x) sapply(cons, function(y) x %in% y))
    # TODO append col_name to each categorical variable name
    colnames(bin)=sapply(f, function(x) paste0(col_name,"_",x))
    j=which(colnames(df) == col_name)
    df2=cbind(df[1:(j-1)], bin*1, df[(j+1):dim(df)[2]])
    return(df2)
}

df=read.table(tsv, fill=TRUE, header=1)
categ=read.table("config/categorical",row.names=1,header=1,fill=TRUE)
for (col in row.names(categ)) {
  drop_cols = unlist(strsplit(categ[col],","))
  df=to_cat(df,col,drop)
}
write.table(df,tsv_out,quote=FALSE,row.names=FALSE)

