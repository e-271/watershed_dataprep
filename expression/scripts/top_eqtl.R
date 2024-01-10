library(dplyr)

# TODO config for the column names of this file.
args = commandArgs(trailingOnly=TRUE)
txt=args[1]

df = read.table(txt, header=T)
df <- df %>% 
        group_by(feature) %>% 
        slice(which.min(pvalue)) %>% 
        select(chr,snp_pos,snp_pos2,ref,alt,feature,pvalue,beta,se,p_het)

write.table(df, "", sep='\t', quote=F, row.names=F, col.names=F)

