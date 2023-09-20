library(dplyr)
library(tidyr)

args = commandArgs(trailingOnly=TRUE)
 
tsv=args[1]
aggf = read.table(args[2], header=F, sep=',', col.names=c("col", "fn"))
type_df = read.table(args[3], header=F, sep=',', row.names=1)

types=as.list(t(type_df))
names(types) = row.names(type_df)

# Read all variants
# For memory efficiency, these can be split by subject & chromosome.
df=read.table(tsv, header=T, sep="\t", colClasses=types)  

# Split values within a column by & or ,
spl <- function(x) {unlist(strsplit(as.character(x),"[,&]"))}
# Standardize NA format
clean <- function(x) { replace(x, x=="." | x=="" | x=="NA", NA) }
# Custom aggregation functions
summ_unique <- function (x) { paste(unique(x[!is.na(x)]), collapse=",")}
summ_append <- function (x) { paste(x[!is.na(x)], collapse=",")}

# Summarise over (subject,gene) pairs using functions defined in aggf
col_max <- aggf %>% filter(fn=="max") %>% select("col")
col_min <- aggf %>% filter(fn=="min") %>% select("col")
col_app <- aggf %>% filter(fn=="append") %>% select("col")
col_uniq <- aggf %>% filter(fn=="unique") %>% select("col")

# Summarise over (subject,gene) pairs
agg_df <- df %>% 
    group_by(Sample,Gene) %>% 
    summarize(across(as.vector(unlist(col_max)), ~max(clean(spl(.x)))), 
              across(as.vector(unlist(col_min)), ~min(clean(spl(.x)))), 
              across(as.vector(unlist(col_app)), ~summ_append(clean(spl(.x)))),
              across(as.vector(unlist(col_uniq)), ~summ_unique(clean(spl(.x)))),
    .groups="keep")

# reorder column names
agg_df <- agg_df[,colnames(df)]

# Write dataframe to stdout
write.table(agg_df,"",quote=FALSE,row.names=FALSE,na="",sep="\t")

