library(tibble)
library(tidyr)
library(dplyr)

args = commandArgs(trailingOnly=TRUE)
eqtls=args[1]

df_eqtls = read.table(eqtls, header=F, col.names=c("SubjectID", "GeneName", "Genotype"))

# Convert genotype to eQTL allele count
# Note as.numeric will induce some NAs for samples with missing genotypes or no eQTL.
df_eqtls <- df_eqtls %>% filter(GeneName != ".") %>%
            separate_wider_position(Genotype,c(GT1=1, 1, GT2=1)) %>%
            mutate(eQTL=as.numeric(GT1) + as.numeric(GT2)) %>%
            select(-GT1, -GT2)
df_eqtls <- df_eqtls %>% pivot_wider(id_cols=GeneName, names_from=SubjectID, values_from=eQTL)
df_eqtls <- df_eqtls %>% column_to_rownames("GeneName")
df_eqtls =  data.frame(t(df_eqtls))

write.table(df_eqtls,"",row.names=T,quote=F)


