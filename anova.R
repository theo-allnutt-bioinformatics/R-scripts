#!/usr/bin/env Rscript

rm(list=ls()) # clear all

args <- commandArgs(trailingOnly=TRUE)
setwd("./")

#anova
data1 <- read.table(args[1],sep="\t", header =TRUE)

meta <- read.table(args[2],sep="\t", header =TRUE)

data1<-cbind(data1,meta)


result <- summary(aov(richness ~ get(args[3]), data=data1))

result <- summary(aov(richness ~ get(args[3]), data=data1))
write("richness",file="alphatests.txt",append=TRUE)
capture.output(result, file="alphatests.txt",append=TRUE)
result <- summary(aov(shannon_2 ~ get(args[3]), data=data1))
write("shannon_2",file="alphatests.txt",append=TRUE)
capture.output(result, file="alphatests.txt",append=TRUE)
result <- summary(aov(simpson ~ get(args[3]), data=data1))
write("simpson",file="alphatests.txt",append=TRUE)
capture.output(result, file="alphatests.txt",append=TRUE)

#k-w
#data <- read.table("~/metagenomics/alpha_grp.tab",sep="\t", header =TRUE)
#result <- summary(kruskal.test(richness ~ group, data=data))
#write("richness",file="~/metagenomics/kw.txt",append=TRUE)
#capture.output(result, file="~/metagenomics/kw.txt",append=TRUE)
