#!/usr/bin/env Rscript
library(dplyr)
library(tidyverse)

args = commandArgs(trailingOnly=TRUE)
tsv=args[1]

df=read.table(tsv, fill=TRUE, header=1,sep="\t",na=c("", "NA"))
pair_df <- df %>% 
            group_by(Gene,POS_ALT) %>% 
            summarise(Sample=paste(Sample,collapse=","), PairN=n(), .groups="keep") %>% 
            filter(PairN > 1) %>% # filter to N>2 
            ungroup() %>% mutate(Pair=row_number()) %>% # number each pair
            separate_longer_delim(Sample, delim=",") %>% # split to 1 row per sample
            reframe(across(Gene:PairN, ~ .x[1:2]), .by=Pair) # remove N>2 pair members (must do this for Watershed to split test set correctly)

df <- df %>% 
        left_join(pair_df) %>% # add pair labels to original dataset
        arrange(Pair) %>% # must sort by pair for Watershed to split test set correctly
        select(-c(POS_ALT)) # drop pos/alt columns

write.table(df,"",quote=FALSE,row.names=FALSE,na="",sep="\t")

