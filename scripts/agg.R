library(dplyr)
library(tidyr)

args = commandArgs(trailingOnly=TRUE)
 
tsv=args[1]
gencode=args[2]
aggf=args[3]
types=args[4]
window=as.numeric(args[5])

# Read all variants
# For memory efficiency, these can be split by subject & chromosome.
type_df = read.table(types, header=F, sep=',', row.names=1)
types=as.list(t(type_df))
names(types) = row.names(type_df)
df=read.table(tsv, header=T, sep="\t", colClasses=types)  
df <- df %>% rename(vepGene = GeneName) # Each variant may have several VEP anno for different genes

aggf = read.table(aggf, header=F, sep=',', col.names=c("col", "fn"))
gencode=read.table(gencode)

# Label gencode gene regions within +-10kb window
colnames(gencode)=c("chr", "start", "end", "strand", "GeneName", "gene_name", "type")
gencode <- gencode %>% 
            select(-strand) %>% 
            mutate(window_start = start-window, window_end = end+window) 
# Inner join filters out variants that are not in a protein-coding / lincRNA gencode gene region
df <- df %>% inner_join(gencode, by=join_by(CHROM==chr, between(POS, window_start, window_end)))
df <- df %>% 
        mutate(distTSS=abs(POS-start), distTES=abs(POS-end)) %>%
        select(-gene_name, -start, -end, -type)

# Reformat pos & alt columns to POS:ALT
df <- df %>% 
        mutate(POS_ALT = paste(POS,ALT,sep=":")) %>% 
        select(-POS,-ALT)

# Custom aggregation functions
# Split values within a column by & or ,
spl <- function(x) {unlist(strsplit(as.character(x),"[,&]"))}
# Standardize NA format
clean <- function(x) { replace(x, x=="." | x=="" | x=="NA", NA) }
# List of comma-separated sorted unique values
summ_unique <- function (x) { paste(sort(unique(x[!is.na(x)])), collapse=",")}

# Summarise over (subject,gene) pairs using functions defined in aggf
col_max <- aggf %>% filter(fn=="max") %>% select("col")
col_min <- aggf %>% filter(fn=="min") %>% select("col")
col_uniq <- aggf %>% filter(fn=="unique") %>% select("col")
col_vep_uniq <- aggf %>% filter(fn=="vep_unique") %>% select("col")

# Summarise over (subject,gene) pairs
agg_df <- df %>% 
    group_by(SubjectID,GeneName) %>% 
    summarize(across(as.vector(unlist(col_max)), ~max(clean(spl(.x)))), 
              across(as.vector(unlist(col_min)), ~min(clean(spl(.x)))), 
              across(as.vector(unlist(col_uniq)), ~summ_unique(clean(spl(.x)))),
              across(as.vector(unlist(col_vep_uniq)), ~summ_unique(clean(spl(.x[GeneName == vepGene] )))),
    .groups="keep")

# reorder column names
agg_df <- agg_df[c("SubjectID", "GeneName", aggf$col)]

# Write dataframe to stdout
write.table(agg_df,"",quote=FALSE,row.names=FALSE,na="",sep="\t")

