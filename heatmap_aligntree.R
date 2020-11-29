#!/usr/bin/Rscript
#source("https://bioconductor.org/biocLite.R")
#biocLite("Heatplus")
#biocLite("vegan")
#biocLite("ape")
#install.packages("phangorn")

#makes a heatmap that has rows and columns ordered by input fasta alignment

library(Heatplus)
library(vegan)
library(RColorBrewer)
library(gplots)
library(ape)
library(phangorn)
library(seqinr)

args = commandArgs(trailingOnly=TRUE)

all.data <- read.table(args[1],sep="\t",header=TRUE,check.names=FALSE)

row.names(all.data) <- all.data$otu
all.data <- all.data[,-1]

head(all.data)

#make data proportions, not abundances, normalise by sample total
#log transfrom if required.
data.prop <- (sweep(all.data, 2, colSums(all.data), FUN="/"))

if (args[4]=="10"){data.s <- log10(data.prop+1)}
if (args[4]=="e"){data.s <- log(data.prop+1)}
if (args[4]=="2"){data.s <- log2(data.prop+1)}
if (args[4]=="a"){data.s <- all.data}
if (args[4]=="p"){data.s <- data.prop}
if (args[4]=="a10"){data.s <- log10(all.data+1)}
if (args[4]=="ae"){data.s <- log(all.data+1)}
if (args[4]=="a2"){data.s <- log2(all.data+1)}

#set group data colours
#group1 <- replace(group1,which(group1==1),"red")
#group1 <- replace(group1,which(group1==2),"blue")
#group1 <-t(group1)
#cbind(names(data.prop), group1)

s1=args[3]
colscale <- colorRampPalette(c("chartreuse4","yellow", "red"), space = "rgb")(100)
x <-data.s[1:s1,]

#make trees
aln <- read.dna(file=args[9],format="fasta")
aln_phy <- phyDat(aln,type="DNA", levels=NULL)
dna_dist <- dist.ml(aln_phy,model="JC")
tree1 <-upgma(dna_dist)
otu.clus <- as.hclust(tree1)

aln2 <- read.dna(file=args[8],format="fasta")
aln_phy2 <- phyDat(aln2,type="DNA", levels=NULL)
dna_dist2 <- dist.ml(aln_phy2,model="JC")
tree2 <-upgma(dna_dist2)
samples.clus <- as.hclust(tree2)

#plot(as.dendrogram(otu.clus), horiz=TRUE, leaflab = "none")

pdf(file=args[2],width=20,height=20)

#set margins
xm = 10 #as.numeric(args[5])
ym = 10 #as.numeric(args[6])

pdf(file=args[2],width=as.numeric(args[5]),height=as.numeric(args[6]))

font_size=as.numeric(args[7])
#mode(xm)
#mode(ym)
heatmap.2(as.matrix(x),Rowv = as.dendrogram(otu.clus), Colv = as.dendrogram(samples.clus), col = colscale,margins=c(xm,ym),trace=c("none"),srtCol=90,key=FALSE, cexRow=font_size, cexCol=font_size)
dev.off()

hm<-heatmap.2(as.matrix(x),Rowv = as.dendrogram(otu.clus), Colv = as.dendrogram(samples.clus), col = colscale,margins=c(xm,ym),trace=c("none"),srtCol=90,key=FALSE)

#order1 <- order.dendrogram(as.dendrogram(otu.clus))

#sorted <- x[match(rev(labels(hm$otu.clus)), rownames(x)), ]

order1 <- x[rev(hm$rowInd), hm$colInd]

write.table(order1,"table-row-order.txt",sep="\t")

