#!/usr/bin/env Rscript

#written for r18 Elcho - endia study


#run metagenome beta diversity analyses
#input is a normalised, e.g. cpm, otu table and a metadata table
#usage:
#pco.R s_.cpm s_


library(vegan)
library("phyloseq")


rm(list=ls()) # clear all

args = commandArgs(trailingOnly=TRUE)
setwd("./")
getwd()

datafile<- args[1]  

outfile<- args[2]

otus <- read.table(file = datafile, header=TRUE, sep="\t",row.names=1,check.names=FALSE,stringsAsFactors = F)


otus<-otu_table(otus,taxa_are_rows=TRUE)

psdata<-phyloseq(otu_table(otus))

pco <- ordinate(psdata, method = "PCoA", distance = "bray")

ev<-pco$values[1]

minev=function(x){x+abs(min(ev))}

evnorm<-sapply(ev,minev)

pc1<-evnorm[1]/sum(evnorm)*100
pc2<-evnorm[2]/sum(evnorm)*100
pc3<-evnorm[3]/sum(evnorm)*100

pcvar<-paste(pc1,pc2,pc3,sep='\t')

print(pcvar)

output1<-pco$vectors[,1:3]

print(output1)

write.table(pcvar,paste(outfile,".pco",sep=""),quote=FALSE,sep="\t")
write.table(output1,paste(outfile,".pco",sep=""),col.names=FALSE,quote=FALSE,sep="\t",append=TRUE)

