#!/usr/bin/env Rscript

#if (!requireNamespace("BiocManager", quietly = TRUE))
    #install.packages("BiocManager")

#BiocManager::install("ALDEx2")

library("ALDEx2")

args = commandArgs(trailingOnly=TRUE)

infile1<-args[1]
meta1<-args[2]
factor1<-args[3]
outfile1<-args[4]

setwd("./")

otu_table <- read.table(file = infile1, header=TRUE, sep="\t",row.names=1,check.names=FALSE,stringsAsFactors = F)

metadata <- read.table(file = meta1, header=TRUE, sep="\t", row.names=1,check.names=FALSE) 


group<-metadata[factor1]
group<-group[,1]

tt <- aldex(otu_table, group,test = "t", mc.samples = 128)

aldex.plot(tt, type = "MW", cutoff = .05 / ncol(otu_table))

write.table(tt,outfile1,sep="\t")
