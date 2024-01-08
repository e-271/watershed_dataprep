library(tidyr)
library(tibble)
library(dplyr)


args = commandArgs(trailingOnly=TRUE)
pv=args[1]
cg=args[2]

pv <- read.table(pv, header=T, sep='\t') 
cg <- read.table(cg, header=F, sep='\t',col.names=c('cluster', 'gene')) 

pv <- pv %>% rename(cluster=CLUSTER_ID)
# Add gene labels to pvalue df
pv <- pv %>% left_join(cg) %>% filter(gene != '.') #df %>% select(cluster,gene) %>% distinct()))

# minimum pvalue
mpv <- pv %>% group_by(gene) %>%  summarize(across(where(is.numeric), min), nc=n())

# Needs to be adjusted later for a given pvalue threshold, using new threshold of 1 - (1 - z)^nc.
write.table(mpv, "", row.names=F, sep='\t', quote=F)




