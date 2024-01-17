library(dplyr)
library(data.table)
library(tibble)
library(tidyr)

args = commandArgs(trailingOnly=TRUE)
cov=args[1]
cts=args[2]
eqtls=args[3]

df_cov <- read.table(cov, header=T, row.names=1)
df_ct <- data.frame(t(read.csv(cts, header=T, row.names=1)))
df_ct = df_ct[row.names(df_cov),]
df_eqtls <- read.table(eqtls, header=T, row.names=1)
df_eqtls = df_eqtls[row.names(df_cov),]
df_eqtls = scale(df_eqtls)

# Check overlap between counts matrix and eQTLs.
cts_genes = colnames(df_ct)
eqtls_genes = colnames(df_eqtls)
shared_genes = intersect(cts_genes,eqtls_genes)
warning(sprintf("%d/%d of RNA-seq measured protein-coding/lincRNA genes have a top eQTL.", length(shared_genes), length(cts_genes)))

# initialize residual matrix
df_resid = data.frame(df_ct)
df_resid[,] = 0

# calculate residuals for each gene
for (gene in colnames(df_ct)) {   
    m = cbind(df_cov, df_ct  %>% select(gene)) %>% rename(geneCounts=gene)
    # residualize eQTL if present for this gene
    e = df_eqtls %>% select(any_of(gene))
    if (dim(e)[2] & !any(is.na(e))) { m = cbind(e %>% rename(eQTL.genotype=gene), m) } 
    model = lm(geneCounts ~ ., data = m)
    df_resid[gene] = model$residuals
    if (any(is.na(model$coefficients))) {
        warning(sprintf("NA coefficients in model residuals for gene %s", gene))
    }
}

# Scale residuals for each gene
df_resid = data.frame(t(scale(df_resid)))

# Reshape to 1 sample to row, and remove Ensembl gene version
df_resid <-  df_resid %>% 
        rownames_to_column("GeneName") %>% 
        separate_wider_delim(GeneName, delim='.', names=c("GeneName", "version")) %>%
        select(-version) %>%
        pivot_longer(!GeneName, names_to="SubjectID", values_to="eOutliers") %>%
        select(SubjectID, GeneName, eOutliers)

# Output residuals
write.table(df_resid,"",row.names=F,quote=F,sep='\t')

