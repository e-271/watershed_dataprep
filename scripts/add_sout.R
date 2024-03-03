library(tidyr)
library(dplyr)

args = commandArgs(trailingOnly=TRUE)
 
tsv=args[1]
sout=args[2]
pvalue=as.numeric(args[3])
nout_std_thresh=as.numeric(args[4])

df=read.table(tsv, header=T, sep="\t", na="")
sdf=read.table(sout, header=T, sep="\t", na="", check.names=F)

# Calulate adjusted pvalue threshold by # clusters/gene
sdf <- sdf %>% 
        mutate(adj_pthresh = pvalue * pvalue / (1 - (1 - pvalue)^nc)) %>% 
        select(-nc)

# Reshape to 1 gene/sample per row
sdf <- sdf %>% 
        pivot_longer(!GeneName & !adj_pthresh, names_to="SubjectID", values_to="sOutliers")

# Adjust sOutlier pvalues by adjusted threshold
sdf <- sdf %>% 
        mutate(adj_sOutliers = sOutliers / (adj_pthresh / pvalue)) %>% 
        select(GeneName, SubjectID, adj_sOutliers) %>% rename(sOutliers=adj_sOutliers)

# Calculate # outliers per sample
df_noutliers <- sdf %>% group_by(SubjectID) %>%
    summarise(noutliers=sum(sOutliers < pvalue))

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

# Add sOutlier adjusted pvalues to tsv (allow missing sOutlier values)
df <- df %>% left_join(sdf)

# Print outlier proportion after filtering
df <- df %>% mutate(sOut = abs(sOutliers) < pvalue)
prop = df %>% summarize(prop = sum(sOut) / n())
df <- df %>% select(-sOut)
warning(sprintf("filtered outlier proportion: %.4f", prop[1,1]))

write.table(df,"",quote=FALSE,row.names=FALSE,na="",sep="\t")

