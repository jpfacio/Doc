#!/bin/bash

seqkit stats -a -T Data/Chen_Analysis/Bins/*.fa.gz \
| awk -v OFS="\t" '
NR==1 {print; next}
{
    sub(/^.*\//,"",$1);      
    sub(/\.fa\.gz$/,"",$1);  
    print
}
' \
> Data/Chen_Analysis/summary.tsv


