#!/stornext/System/data/apps/R/R-3.5.1/lib64/R/bin/Rscript

#source("https://bioconductor.org/biocLite.R")
#biocLite("Heatplus")
#biocLite("vegan")
#biocLite("ape")
#install.packages("phangorn")
#install.packages("phytools")

#makes a heatmap that has rows ordered by input tree

library(Heatplus)
library(vegan)
library(RColorBrewer)
library(gplots)
library(ape)
library(phangorn)
library(seqinr)
library(phytools)

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
#aln <- read.dna(file=args[9],format="fasta")
#aln_phy <- phyDat(aln,type="DNA", levels=NULL)
#dna_dist <- dist.ml(aln_phy,model="JC")

tree1 <-read.tree(file=args[5])

tree1<- chronos(tree1)
tree1<-root(tree1)
tree1=midpoint.root(tree1)

print (tree1)
is.ultrametric(tree1)
is.rooted(tree1)

otu.clus <- as.hclust(tree1)


#aln2 <- read.dna(file=args[8],format="fasta")
#aln_phy2 <- phyDat(aln2,type="DNA", levels=NULL)
#dna_dist2 <- dist.ml(aln_phy2,model="JC")
#tree2 <-upgma(dna_dist2)
#samples.clus <- as.hclust(tree2)

#plot(as.dendrogram(otu.clus), horiz=TRUE, leaflab = "none")

pdf(file=args[2],width=20,height=20)

#set margins
xm = 10 #as.numeric(args[6])
ym = 10 #as.numeric(args[7])

font_size=as.numeric(args[8])
#mode(xm)
#mode(ym)
heatmap.2(as.matrix(x),Rowv = as.dendrogram(otu.clus), Colv=FALSE,col = colscale,margins=c(xm,ym),trace=c("none"),srtCol=90,key=FALSE, cexRow=font_size, cexCol=font_size)
dev.off()

hm<-heatmap.2(as.matrix(x),Rowv = as.dendrogram(otu.clus), col = colscale,margins=c(xm,ym),trace=c("none"),srtCol=90,key=FALSE)

#order1 <- order.dendrogram(as.dendrogram(otu.clus))

#sorted <- x[match(rev(labels(hm$otu.clus)), rownames(x)), ]

order1 <- x[rev(hm$rowInd), hm$colInd]

write.table(order1,"table-row-order.txt",sep="\t")

