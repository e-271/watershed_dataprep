#!/usr/bin/env Rscript
library(dplyr)
library(tidyverse)

args = commandArgs(trailingOnly=TRUE)
tsv=args[1]
pairs=args[2]

df=read.table(tsv, fill=TRUE, header=1,sep="\t",na=c("", "NA"))
pair_df=read.table(pairs, fill=TRUE, header=1,sep="\t",na=c("", "NA")) %>% select(-AF) 

# Add pos/alt from pair df (includes indels / infrequent variants)
pair_df <- pair_df %>% rename(POS_ALT_pair = POS_ALT)
df <- df %>% left_join(pair_df)

# Construct pair dataset
pairs <- df %>% 
            group_by(GeneName,POS_ALT_pair) %>%
            summarise(SubjectID=paste(SubjectID,collapse=","), PairN=n(), .groups="keep") %>%
            filter(PairN > 1) %>% # filter to N>2
            ungroup() %>% mutate(N2pair=row_number()) %>% # number each pair
            separate_longer_delim(SubjectID, delim=",") %>% # split to 1 row per sample
            reframe(across(GeneName:PairN, ~ .x), .by=N2pair)
warning(sprintf("found %d gene/subjects in pairs", dim(pairs)[1]))

# Label df by pairs and order by pair number
df <- df %>%
        left_join(pairs) %>%
        arrange(N2pair) %>% 
        select(-POS_ALT)

# Make training set from non-pairs / pairs that became n=1 after filtering
train = df %>% filter(is.na(N2pair)) #%>% mutate(N2pair=na_if(PairN, 1)) %>% select(-PairN)

# Make test set from N>=2 pairs
# Drop N>3 samples but keep the first 2
#test_nplus = df %>% filter(!is.na(N2pair)) %>%
#                group_by(N2pair, GeneName) %>% 
#                summarize(N1=SubjectID[1], N2=SubjectID[2])
#stopifnot(sum(is.na(test_nplus$N2)) == 0) # make sure every pair has an N2 member
#test_nplus <- test_nplus %>%
#                pivot_longer(c("N1","N2"), values_to="SubjectID") %>% 
#                select(-name)
#test = df %>% inner_join(test_nplus)
# Make test set from N=2 pairs only
test = df %>% filter(PairN == 2)

df <- rbind(train, test) %>% select(-POS_ALT_pair, -PairN)

write.table(df,"",quote=FALSE,row.names=FALSE,na="",sep="\t")

