#!/usr/bin/env Rscript
library(dplyr)
library(tidyverse)

args = commandArgs(trailingOnly=TRUE)
tsv=args[1]

df=read.table(tsv, fill=TRUE, header=1,sep="\t",na=c("", "NA"))
pair_df <- df %>% 
            group_by(GeneName,POS_ALT) %>% 
            summarise(SubjectID=paste(SubjectID,collapse=","), PairN=n(), .groups="keep") %>% 
            filter(PairN > 1) %>% # filter to N>2 
            ungroup() %>% mutate(N2pair=row_number()) %>% # number each pair
            separate_longer_delim(SubjectID, delim=",") %>% # split to 1 row per sample
            reframe(across(GeneName:PairN, ~ .x[1:2]), .by=N2pair) # remove N>2 pair members (must do this for Watershed to split test set correctly)

df <- df %>% 
        left_join(pair_df) %>% # add pair labels to original dataset
        arrange(N2pair) %>% # must sort by pair for Watershed to split test set correctly
        select(-POS_ALT, -PairN) # drop pos/alt and PairN columns

write.table(df,"",quote=FALSE,row.names=FALSE,na="",sep="\t")

