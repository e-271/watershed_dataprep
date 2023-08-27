#!/bin/bash

alias max="sort -n | tail -n 1"
alias min="sort -n | head -n 1"
alias unique="sort -u | tr '\n' ',' | sed 's/,$//'"
alias append="tr '\n' ',' | sed 's/,$//'"
alias split="sed 's/[,&]/\n/g'"
alias clean_missing="sed 's/^\.$//g' | sed 's/^NA$//g'"

