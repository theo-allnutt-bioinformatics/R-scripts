#!/usr/bin/Rscript
#biocLite("Heatplus")
#biocLite("vegan")
#biocLite("ape")

#makes a heatmap that has rows ordered by tree of input alignment

library(Heatplus)
library(vegan)
library(RColorBrewer)
library(gplots)
library(ape)
library(phangorn)
library(seqinr)

args = commandArgs(trailingOnly=TRUE)

all.data <- read.table(args[1],sep="\t",header=TRUE)

row.names(all.data) <- all.data$otu
all.data <- all.data[,-1]
group1 = all.data[1,]
all.data <- all.data[-1,]
head(all.data)

#make data proportions, not abundances, normalise by sample total
#log transfrom if required.
data.prop <- (sweep(all.data, 2, colSums(all.data), FUN="/"))

if (args[5]=="10"){data.s <- log10(data.prop+1)}
if (args[5]=="e"){data.s <- log(data.prop+1)}
if (args[5]=="2"){data.s <- log2(data.prop+1)}
if (args[5]=="a"){data.s <- all.data}
if (args[5]=="p"){data.s <- data.prop}
if (args[5]=="a10"){data.s <- log10(all.data+1)}
if (args[5]=="ae"){data.s <- log(all.data+1)}
if (args[5]=="a2"){data.s <- log2(all.data+1)}


#set group data colours
group1 <- replace(group1,which(group1==1),"black")
group1 <- replace(group1,which(group1==2),"grey") #add a  line here to add groups and colours
group1 <- replace(group1,which(group1==3),"blue")
group1 <- replace(group1,which(group1==4),"red")
group1 <- replace(group1,which(group1==5),"yellow")
group1 <- replace(group1,which(group1==6),"orange")
group1 <- replace(group1,which(group1==7),"purple")
group1 <- replace(group1,which(group1==8),"green")
group1 <-t(group1)
#cbind(names(data.prop), group1)

s1=args[3]
colscale <- colorRampPalette(c("chartreuse4","yellow", "red"), space = "rgb")(100)

#simple map
#heatmap(as.matrix(data.prop[1:s1,]), Rowv = NA, Colv = NA, col = scaleyellowred, margins = c(10, 2))

#make trees
x <-data.s[1:s1,]
#otu.dist <- vegdist(x, method = "bray")
#tree1 <- read.tree(file=args[4])#
#utree = chronos(tree1, lambda = 0, model = "correlated")
#otu.clus <- as.hclust(utree)
aln <- read.dna(file=args[4],format="fasta")
aln_phy <- phyDat(aln,type="DNA", levels=NULL)

dna_dist <- dist.ml(aln_phy,model="JC")
tree1 <-upgma(dna_dist)
otu.clus <- as.hclust(tree1)

samples.dist <- vegdist(t(data.prop), method = "bray")
samples.clus <- hclust(samples.dist, "aver")

#plot(as.dendrogram(otu.clus), horiz=TRUE)
pdf(file=args[2],width=10,height=12)

heatmap.2(as.matrix(x),Rowv = as.dendrogram(otu.clus), Colv = as.dendrogram(samples.clus), col = colscale,margins=c(3,12), ColSideColors = group1,trace=c("none"),srtCol=90,key=FALSE)
dev.off()

heatmap.2(as.matrix(x),Rowv = as.dendrogram(otu.clus), Colv = as.dendrogram(samples.clus), col = colscale,margins=c(3,12), ColSideColors = group1,trace=c("none"),srtCol=90,key=FALSE)

hm<-heatmap.2(as.matrix(x),Rowv = as.dendrogram(otu.clus), Colv = as.dendrogram(samples.clus), col = colscale,margins=c(3,12), ColSideColors = group1,trace=c("none"),srtCol=0,key=FALSE)

sorted <- x[match(rev(labels(hm$rowDendrogram)), rownames(x)), ]

write.table(sorted,"table.txt")

pdf(file="plot.pdf",width=7,height=50)
plot(tree1,font=1,direction="rightwards",cex=0.1,show.tip.label = FALSE)
dev.off()
s2<-as.numeric(s1)
h1<-s2*50
jpeg(file="plot.jpg",width=720,height=h1)
plot(tree1,font=1,direction="rightwards",cex=0.1,show.tip.label = FALSE, no.margin=TRUE,edge.width = 5)
dev.off()
