library(dplyr)
library(tibble)

args = commandArgs(trailingOnly=TRUE)
meta=args[1]
cfg=args[2]

df=read.csv(meta,header=T)
cdf=read.csv(cfg,header=T)

# Subset to covariates in the config file
df <- df %>% select(cdf$Name)

# Set sample column to DF row name.
sample_id = deframe(cdf %>% filter(Type=='Sample') %>% select(Name))
df <- df %>% column_to_rownames(var = sample_id)

# Get list of covariates of each type
bin_covs = deframe(cdf %>% filter(Type=='Binary') %>% select(Name))
cat_covs = deframe(cdf %>% filter(Type=='Categorical') %>% select(Name))

# Cast covariates to Binary (convert to Factor then to Numeric)
df <- df %>% mutate(across(bin_covs, ~as.numeric(as.factor(.x))-1))
if (!all(df %>% select(bin_covs) < 2)) {warning("WARNING: Found more than 2 factors for a binary variable.")}

# Cast covariates to Categorical
to_categorical <- function(df, col) {
    # categorical conversion. drops the 1st batch label to preserve matrix rank 
    nl = length(levels(df[,col]))
    cat <- sapply(levels(df[,col])[2:nl], function(x) {as.numeric(x == df[,col])} ) 
    # rename categories with 'col.' prefix             
    colnames(cat) <- paste0(col,".",colnames(cat)) # rename
    # return original df with 'col' replaced by categorical
    return(cbind(df %>% select(-col), cat))
}
df <- df %>% mutate(across(cat_covs, ~as.factor(.x)))
cats <- lapply(df %>% select(cat_covs), levels)
for (cat in cat_covs) {
    df <- to_categorical(df, cat)   
}

# Filter out linearly dependent covariates.
# source: https://stackoverflow.com/questions/19100600/extract-maximal-set-of-independent-columns-from-a-matrix
d1 = dim(df)[2]
df <-  df[, qr(df)$pivot[seq_len(qr(df)$rank)]]
d2 = dim(df)[2]
warning(sprintf("Removed %d/%d linearly dependent columns.", d1-d2, d1))

# Output covariates
write.table(df,"",row.names=T,quote=F,sep='\t')



