library(tidyr)
library(dplyr)

args = commandArgs(trailingOnly=TRUE)
 
tsv=args[1]
sout=args[2]
pvalue=as.numeric(args[3])

df=read.table(tsv, header=T, sep="\t", na="")
sdf=read.table(sout, header=T, sep="\t", na="")

# Calulate adjusted pvalue threshold by # clusters/gene
sdf <- sdf %>% 
        mutate(adj_pthresh = pvalue * pvalue / (1 - (1 - pvalue)^nc)) %>% 
        select(-nc)

# Reshape to 1 gene/sample per row
sdf <- sdf %>% 
        pivot_longer(!gene & !adj_pthresh, names_to="SubjectID", values_to="sOutliers") %>% 
        rename(GeneName=gene)

# Adjust sOutlier pvalues by adjusted threshold
sdf <- sdf %>% 
        mutate(adj_sOutliers = sOutliers / (adj_pthresh / pvalue)) %>% 
        select(GeneName, SubjectID, adj_sOutliers) %>% rename(sOutliers=adj_sOutliers)

# Add sOutlier adjusted pvalues to tsv (allow missing sOutlier values)
df <- df %>% left_join(sdf)

write.table(df,"",quote=FALSE,row.names=FALSE,na="",sep="\t")

