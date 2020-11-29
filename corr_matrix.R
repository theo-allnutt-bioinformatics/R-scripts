#!/usr/bin/Rscript


#makes a correlation matrix from two otu input tables

library(Hmisc)
library(vcd)

args = commandArgs(trailingOnly=TRUE)

x <- read.table(args[1],sep="\t",header=TRUE,check.names=FALSE)

y <- read.table(args[2],sep="\t",header=TRUE,check.names=FALSE)

a<-unlist(strsplit(args[1],"\\."))[1]

b<-unlist(strsplit(args[2],"\\."))[1]

name1 <- paste(a,"_",b,"_pearson.tab",sep="")
name2 <- paste(a,"_",b,"_spearman.tab",sep="")

row.names(x) <- x$otus
x <- x[,-1]
row.names(y) <- y$otus
y <- y[,-1]

x <- t(x)

y <- t(y)

head(x)

head(y)

c1 <- rcorr(x,y, type="pearson")

c2 <- rcorr(x,y, type="spearman")

c1
c2

c1.r <- data.frame(c1$r)
c1.p <- data.frame(c1$P)

c2.r <- data.frame(c2$r)
c2.p <- data.frame(c2$P)

write.table("Pearson R\n",name1,sep="\t",quote=FALSE,col.names=TRUE,row.names=TRUE)
write.table(c1.r,name1,append=TRUE,sep="\t",quote=FALSE,col.names=TRUE,row.names=TRUE)
write.table("Pearson P-value\n",name1,append=TRUE,sep="\t",quote=FALSE,col.names=TRUE,row.names=TRUE)
write.table(c1.p,name1,append=TRUE,sep="\t",quote=FALSE,col.names=TRUE,row.names=TRUE)

write.table("Spearman R\n",name2,sep="\t",quote=FALSE,col.names=TRUE,row.names=TRUE)
write.table(c2.r,name2,append=TRUE,sep="\t",quote=FALSE,col.names=TRUE,row.names=TRUE)
write.table("Spearman P-value\n",name2,append=TRUE,sep="\t",quote=FALSE,col.names=TRUE,row.names=TRUE)
write.table(c2.p,name2,append=TRUE,sep="\t",quote=FALSE,col.names=TRUE,row.names=TRUE)

 



