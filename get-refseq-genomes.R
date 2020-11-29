#!/usr/bin/env Rscript



library(biomartr)

args <- commandArgs(trailingOnly=TRUE)
setwd("./")

#make args a variable= get(args[1])

print (args[1])

meta.retrieval(kingdom = args[1], db = "refseq", type = "genome")

