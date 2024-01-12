library(PCAtools)

args = commandArgs(trailingOnly=TRUE)
counts=args[1]

cts <- read.csv(counts, header=T, sep=",", row.names=1)
# Note in the prev. step we scale/center each *gene*, here we scale/center each *sample* for PCA.
p <- pca(cts, scale=T, center=T) 
pcs <- p$rotated

# Output PCs
write.table(pcs,"",row.names=T,quote=F,sep='\t')
