library(dplyr)
library(tidyr)

args = commandArgs(trailingOnly=TRUE)
 
tsv=args[1]
gencode=args[2]
aggf=args[3]
types=args[4]
window=as.numeric(args[5])
match_vep_gene=as.logical(args[6])

# Read column types
type_df = read.table(types, header=F, sep=',', row.names=1)
types=as.list(t(type_df))
names(types) = row.names(type_df)

# Aggregation functions
aggf = read.table(aggf, header=F, sep=',', col.names=c("col", "fn"))
col_max <- aggf %>% filter(fn=="max") %>% select("col") %>% pull()
col_min <- aggf %>% filter(fn=="min") %>% select("col") %>% pull()
col_uniq <- aggf %>% filter(fn=="unique") %>% select("col") %>% pull()

# TSV variants
df=read.table(tsv, header=T, sep="\t", colClasses=types)  
if (match_vep_gene) {
    df <- df %>% rename(vepGene = GeneName) 
}
# Label gencode gene regions within +- window
gencode=read.table(gencode)
colnames(gencode)=c("chr", "start", "end", "strand", "GeneName", "gene_name", "type")
gencode <- gencode %>% 
            select(-strand) %>% 
            mutate(window_start = start-window, window_end = end+window) 

# Inner join filters out variants that are not in a protein-coding / lincRNA gencode gene region
# vepGene == '.' occurs at VEP regulatory regions
df <- df %>% inner_join(gencode, by=join_by(CHROM==chr, between(POS, window_start, window_end)))
if (match_vep_gene) { df <- df %>% filter(vepGene == GeneName | vepGene == ".") }
df <- df %>%
    mutate(distTSS=abs(POS-start), distTES=abs(POS-end)) %>%
    mutate(POS_ALT = paste(POS,ALT,sep=":")) %>% 
    select(-gene_name, -start, -end, -type, -POS, -ALT)

# Split values within a column by & or ,
spl <- function(x) {unlist(strsplit(as.character(x),"[,&]"))}
# Standardize NA format
clean <- function(x) { replace(x, x=="." | x=="" | x=="NA", NA) }
# List of comma-separated sorted unique values
summ_unique <- function (x) { paste(sort(unique(x[!is.na(x)])), collapse=",")}

# Summarise over (subject,gene) pairs
agg_df <- df %>% 
    group_by(SubjectID,GeneName) %>% 
    summarize(across(col_max, ~max(as.numeric(clean(spl(.x))))), 
              across(col_min, ~min(as.numeric(clean(spl(.x))))), 
              across(col_uniq, ~summ_unique(clean(spl(.x)))),
    .groups="keep")

# reorder column names
agg_df <- agg_df[c("SubjectID", "GeneName", aggf$col)]

# Write dataframe to stdout
write.table(agg_df,"",quote=FALSE,row.names=FALSE,na="",sep="\t")

