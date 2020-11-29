#!/usr/bin/Rscript
#biocLite("Heatplus") #uncomment if you need to install these libraries
#biocLite("vegan")
#biocLite("ape")

#Theo Allnutt, 2016.
#makes a heatmap that has columns ordered by upgma of variable with rows in original order - no sorting
#usage:
#$heatmap3.R otu-table.txt output-image.pdf 9 a 3 33
#n.b. otu table must be in text format and have the top row '#' removed from '#otu id'.
#'9'=number of otus to draw starting at top of table
#'a'=raw abundance data from table:
#options for transformation:
#'10'=proportion log10; 'e'=proportion natural log; '2'=proportion log base 2; 
#'p'=proportion=columns divided by sum of column; 'a10'=raw data log10; 
#'ae'=raw data log10; 'a2'=raw data log2.
#'3'=bottom margin for labels
#'33'=right margin for labels

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

if (args[4]=="10"){data.s <- log10(data.prop+1)}
if (args[4]=="e"){data.s <- log(data.prop+1)}
if (args[4]=="2"){data.s <- log2(data.prop+1)}
if (args[4]=="a"){data.s <- all.data}
if (args[4]=="p"){data.s <- data.prop}
if (args[4]=="a10"){data.s <- log10(all.data+1)}
if (args[4]=="ae"){data.s <- log(all.data+1)}
if (args[4]=="a2"){data.s <- log2(all.data+1)}

#set group data colours
group1 <- replace(group1,which(group1==1),"black")
group1 <- replace(group1,which(group1==2),"grey")
group1 <-t(group1)
#cbind(names(data.prop), group1)

s1=args[3]
colscale <- colorRampPalette(c("chartreuse4","yellow", "red"), space = "rgb")(100)

#simple map
#heatmap(as.matrix(data.prop[1:s1,]), Rowv = NA, Colv = NA, col = scaleyellowred, margins = c(10, 2))

#make trees
x <-data.prop[1:s1,]
#otu.dist <- vegdist(x, method = "bray")
#otu.clus <- hclust(otu.dist, "aver")
samples.dist <- vegdist(t(data.prop), method = "bray")
samples.clus <- hclust(samples.dist, "aver")

#plot(as.dendrogram(otu.clus), horiz=TRUE)
pdf(file=args[2],width=12,height=12)
xm = as.numeric(args[5])
ym = as.numeric(args[6])
#mode(xm)
#mode(ym)


lwid = c(0.05,1) #set width left margin
lhei = c(2,as.numeric(s1)*0.75) #set dendrogram, heatmap height

heatmap.2(as.matrix(x),cexRow=2,lwid=lwid,lhei=lhei,dendrogram = c("column"),Rowv=FALSE,Colv = as.dendrogram(samples.clus), col = colscale,margins=c(xm,ym), ColSideColors = group1,trace=c("none"),srtCol=0,key=FALSE)
dev.off()

