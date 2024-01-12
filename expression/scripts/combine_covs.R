library(dplyr)
library(tibble)
library(PCAForQTL)

args = commandArgs(trailingOnly=TRUE)
cov=args[1]
gpc=args[2]
pca=args[3]
gpc_k=as.numeric(args[4])
pc_k=as.numeric(args[5])

cov_df=read.table(cov,header=T,row.names=1) %>% tibble::rownames_to_column("Sample")
gpc_df=read.table(gpc,header=F,row.names=1) %>% tibble::rownames_to_column("Sample")
pca_df=read.table(pca,header=T,row.names=1) %>% tibble::rownames_to_column("Sample")

# Drop unused PCs.
gpc_df <- gpc_df[,1:gpc_k]
pca_df <- pca_df[,1:pc_k]

# Merge hidden covariates
d1 = dim(pca_df)[1]
hidden_df <- pca_df %>% inner_join(gpc_df[,1:gpc_k])
d2 = dim(hidden_df)[1]
hidden_cols = colnames(hidden_df %>% select(-Sample))
warning(sprintf("Filtered out %d/%d samples with missing genetic PCs.", d1-d2, d1))

# Merge hidden & known covariates.
d1 = dim(cov_df)[1]
df <- cov_df %>% inner_join(hidden_df)
d2 = dim(df)[1]
known_cols = colnames(cov_df %>% select(-Sample))
warning(sprintf("Filtered out %d/%d samples with missing hidden covariates.", d1-d2, d1))


# Filter out hidden covariates that have high covariance with known covariates.
df <- df %>% tibble::column_to_rownames("Sample")
known_filt <- PCAForQTL::filterKnownCovariates(df %>% select(all_of(known_cols)),
        df %>% select(all_of(hidden_cols)),
        unadjustedR2_cutoff=0.7, verbose=FALSE)
warning(sprintf("Removed known covariates with high hidden covariate correlation: %s", 
                paste(setdiff(known_cols, colnames(known_filt)))))

# Combine filtered known covs with hidden, & scale each column.
df <- cbind(known_filt, df %>% select(all_of(hidden_cols)))
df <- scale(df)

# Output covariates
write.table(df,"",row.names=T,quote=F,sep='\t')



