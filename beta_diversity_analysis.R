#!/usr/bin/env Rscript

#written for r18 Elcho - endia study

#devtools::install_github("vmikk/metagMisc")

#run metagenome beta diversity analyses
#input is a normalised, e.g. cpm, otu table and a metadata table
#usage:
#beta_diversity_analysis.py s.cpm meta.txt Study
#where "Study" is the fixed variable



library(vegan)
library("phyloseq")
library("pairwiseAdonis")

rm(list=ls()) # clear all

args = commandArgs(trailingOnly=TRUE)
setwd("./")
getwd()

datafile<- args[1] 
metafile<- args[2] 

outfile<- args[3]

otus <- read.table(file = datafile, header=TRUE, sep="\t",row.names=1,check.names=FALSE,stringsAsFactors = F)

metadata <- read.table(file = metafile, header=TRUE, sep="\t", row.names=1,check.names=FALSE) 

var1 <-args[4] 

otus<-otu_table(otus,taxa_are_rows=TRUE)

psdata<-phyloseq(otu_table(otus),sample_data(metadata))

groupvar<-metadata[[var1]]

bc<-distance(psdata,method="bray")

ad<-adonis2(bc ~ groupvar,metadata,method="bray",permutations=999)

print(ad)

filename<-paste(outfile,".adonis")

write.table(ad,file=filename,quote=FALSE,sep="\t")

pw<-pairwise.adonis(bc,groupvar)

print(pw)

write.table(pw,file=filename,col.names=FALSE,quote=FALSE,sep="\t",append=TRUE)




