library(tidyr)
library(dplyr)

args = commandArgs(trailingOnly=TRUE)
 
tsv=args[1]
zthresh=as.numeric(args[2])
nout_std_thresh=as.numeric(args[3])

df=read.table(tsv, header=T, sep="\t", na="")

# Genes with at least 1 outlier
has_outliers <- df %>% group_by(Gene) %>% 
    filter(abs(eOutliers) > zthresh) %>% 
    summarise()
has_outliers <- unlist(has_outliers)

# Calculate # outliers per sample
df_noutliers <- df %>% group_by(Sample) %>%
    summarise(noutliers=sum(abs(eOutliers) > zthresh))

# Calculate # outlier threshold as (mean + std * nout_std_thresh)
nu = mean(df_noutliers$noutliers)
ns = sqrt(var(df_noutliers$noutliers))
nout_thresh = round(nu + nout_std_thresh * ns)
warning(sprintf("using maximum outlier threshold %d", nout_thresh))

# Filter
under_outlier_thresh <- df_noutliers %>%
    filter(noutliers < nout_thresh) %>%
    select(Sample)
under_outlier_thresh <- unlist(under_outlier_thresh)

df_keep <- df %>% filter(Gene %in% has_outliers) %>% filter(Sample %in% under_outlier_thresh)

# Convert zscores to signed pvalues
df_keep = df_keep %>% mutate(eOutliers=sign(eOutliers) * pnorm(-abs(eOutliers))*2)

# Convert zscores to signed pvalues
write.table(df_keep,"",quote=FALSE,row.names=FALSE,na="",sep="\t")

