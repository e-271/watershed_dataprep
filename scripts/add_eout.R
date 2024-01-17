library(tidyr)
library(dplyr)

args = commandArgs(trailingOnly=TRUE)
 
tsv=args[1]
eout=args[2]
pvalue=as.numeric(args[3])
nout_std_thresh=as.numeric(args[4])

zthresh=abs(qnorm(pvalue/2))
df=read.table(tsv, header=T, sep="\t", na="")
edf=read.table(eout, header=T, sep="\t", na="")

# Add eOutlier zscores
df <- df %>% inner_join(edf)

# Filter out genes with no outliers
has_outliers <- df %>% group_by(GeneName) %>% 
    filter(abs(eOutliers) > zthresh) %>% 
    summarise()
has_outliers <- unlist(has_outliers)

d1 = dim(df)[1]
df <- df %>% filter(GeneName %in% has_outliers)
d2 = dim(df)[1]
warning(sprintf("filtered out %d/%d entries by genes with no outliers",d1-d2,d1 ))

# Calculate # outliers per sample
df_noutliers <- df %>% group_by(SubjectID) %>%
    summarise(noutliers=sum(abs(eOutliers) > zthresh))

# Calculate # outlier threshold as (mean + std * nout_std_thresh)
nu = mean(df_noutliers$noutliers)
ns = sqrt(var(df_noutliers$noutliers))
nout_thresh = round(nu + nout_std_thresh * ns)
warning(sprintf("using maximum outlier threshold %d", nout_thresh))

# Filter out individuals over outlier threshold
under_outlier_thresh <- df_noutliers %>%
    filter(noutliers < nout_thresh) %>%
    select(SubjectID)
under_outlier_thresh <- unlist(under_outlier_thresh)

d1 = dim(df)[1]
df <- df %>% filter(SubjectID %in% under_outlier_thresh)
d2 = dim(df)[1]
warning(sprintf("filtered out %d/%d entries by max subject #outlier threshold",d1-d2,d1 ))

# Print outlier proportion after filtering
df <- df %>% mutate(eOut = abs(eOutliers) > zthresh)
prop = df %>% summarize(prop = sum(eOut) / n())
warning(sprintf("filtered outlier proportion: %.4f", prop[1,1]))

# Rescale zscores
# df <- df %>% mutate(eOutliers=scale(eOutliers))

# Convert zscores to signed pvalues
df = df %>% mutate(eOutliers=sign(eOutliers) * pnorm(-abs(eOutliers))*2)

# Print filtered/rescaled outlier proportion
df <- df %>% mutate(eOut = abs(eOutliers) < 0.01)
prop = df %>% summarize(prop = sum(eOut) / n())
df <- df %>% select(-eOut)
warning(sprintf("filtered outlier proportion after rescaling: %.4f", prop[1,1]))

write.table(df,"",quote=FALSE,row.names=FALSE,na="",sep="\t")

