library(dplyr)
library(data.table)

args = commandArgs(trailingOnly=TRUE)
cov=args[1]
cts=args[2]

df_cov <- read.table(cov, header=T, row.names=1)
df_ct <- data.frame(t(read.csv(cts, header=T, row.names=1)))
df_ct = df_ct[row.names(df_cov),]

# initialize residual matrix
df_resid = data.frame(df_ct)
df_resid[,] = 0

# calculate residuals for each gene
for (gene in colnames(df_ct)) {   
    m = cbind(df_cov, df_ct  %>% select(gene)) %>% 
        rename(geneCounts=gene)
    model = lm(geneCounts ~ ., data = m)
    df_resid[gene] = model$residuals
    if (any(is.na(model$coefficients))) {
        warning(sprintf("NA coefficients in model residuals for gene %s", gene))
    }
}

# Scale residuals for each gene
df_resid = data.frame(t(scale(df_resid)))

# Output residuals
write.table(df_resid,"",row.names=T,quote=F,sep='\t')

