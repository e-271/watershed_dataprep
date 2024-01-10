library(edgeR)
library(dplyr)

args = commandArgs(trailingOnly=TRUE)
counts=args[1]
gene_lengths=args[2]
length_type=args[3]
min_ct=as.numeric(args[4])
min_prop=as.numeric(args[5])

cts <- read.csv(counts, header=T, sep=",", row.names=1)
gl <- read.table(gene_lengths, header=T,row.names=1)
gl <- gl %>% select(length_type)
gl <- gl[row.names(cts),]

# Calculate TPM
# following https://www.biostars.org/p/388584/
dge <- DGEList(counts=cts)
dge <- calcNormFactors(dge)
RPKM <- rpkm(dge,log=F, gene.length=gl)
TPM <- t( t(RPKM) / colSums(RPKM) ) * 1e6

# Filter genes where 20% of subjects have TPM <0.1
d1=dim(TPM)[1]
over_proportion <- rowSums(TPM > min_ct) / dim(TPM)[2]
TPM <- TPM[over_proportion > min_prop,]
d2=dim(TPM)[1]
warning(sprintf("filtered out %d/%d genes where >%d%% subjects have <%.1f TPM", d1-d2, d1, round(min_prop*100), min_ct))

# Log-transform the TPM
logtpm = log2(TPM+2)

# Scale expression of each gene to mean 0 variance 1
logtpm.scaled = t(scale(t(logtpm), center=T, scale=T))
cts = rbind(cts$Gene, logtpm.scaled)

# Output normalized counts
write.csv(cts,"",row.names=T,quote=F)

