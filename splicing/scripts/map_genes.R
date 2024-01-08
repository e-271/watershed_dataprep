library(dplyr)
library(tibble)

args = commandArgs(trailingOnly=TRUE)
tsv=args[1]
genes=args[2]

dfg=read.table(genes, sep='\t', header=F, col.names=c("cluster", "gene"))
dfg <- dfg %>% group_by(cluster) %>% summarize(gene = paste(gene,collapse=','))

df = read.table(tsv, sep=' ', row.names=1, header=T)
df["cluster"] = unlist(lapply(row.names(df), function(s) unlist(strsplit(s, ':'))[4] ))

df_gene <- df %>% 
           rownames_to_column(var="rowname") %>%
           inner_join(dfg) %>%
           mutate(rowname=paste(rowname,gene,sep=":")) %>% 
           column_to_rownames(var="rowname") %>% 
           select(-gene,-cluster)

write.table(df_gene, "", sep=' ', quote=F, row.names=T)


