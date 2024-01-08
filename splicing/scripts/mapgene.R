library(dplyr)

tsv= "filtered.AFRIntrons_perind_numers.counts"
df = read.table(tsv, sep=' ', row.names=1, header=T)

# Ben: Removed all exon-exon junctions belonging to a LeafCutter cluster where more than
# 10% of the samples had less than 3 reads summed across all exon-exon junctions
# assigned to that LeafCutter cluster.
df_s = df %>% group_by(cluster) %>% summarize_all(sum)
pct = rowSums((df_s %>% select_if(is.numeric)) < 3) / (dim(df)[2]-1)
keep_cl = pct < 0.1
df_keep <- df %>% filter(cluster %in% unlist(df_s[keep_cl, "cluster"])) %>% select(-cluster)
write.table(df_keep, "cluster_filtered.AFRIntrons_perind_numers.counts", sep=' ', quote=F, row.names=T)

