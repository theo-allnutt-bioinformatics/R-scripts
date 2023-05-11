#!/usr/bin/env Rscript

#usage:
#tree_compare.R tree1 tree2 tree-out.pdf outgroup_name

library(phytools)
library(ggtree)

rm(list=ls()) # clear all

args = commandArgs(trailingOnly=TRUE)

treefile1<-args[1]
treefile2<-args[2]
treeout<-args[3]
outg<-args[4]

wd<-"./"
setwd(wd)
getwd()

print(treefile1)
print(treefile2)

tree1<-ape::read.tree(file=treefile1)

tree2<-ape::read.tree(file=treefile2)

#get outgroup node number

#outg<-'377952_Paraserianthes_lophantha'

c<-0
for (x in tree1$tip.label){
  c=c+1
  if (x==outg){rootnum1<-c}
}

print(rootnum1)

tree1<-reroot(tree1, rootnum1,position=0.5*tree1$edge.length[which(tree1$edge[,2]==rootnum1)])

rootedge.to.singleton(tree1)

c<-0
for (x in tree2$tip.label){
  c=c+1
  if (x==outg){rootnum2<-c}
}

print(rootnum2)

tree2<-reroot(tree2, rootnum2, position=0.5*tree2$edge.length[which(tree2$edge[,2]==rootnum2)])

rootedge.to.singleton(tree2)

ladderize(tree1)

ladderize(tree2)

hp.cophylo<-cophylo(tree1,tree2,rotate=FALSE)  #angtree on left

pdf(file=treeout,width=50,height=80)

plot(hp.cophylo,link.type="curved",link.lwd=2,link.lty="solid",link.col=make.transparent("red",0.6),fsize=2.2)

edgelabels.cophylo(hp.cophylo$trees[[1]]$node.label[2:hp.cophylo$trees[[1]]$Nnode],edge=sapply(2:hp.cophylo$trees[[1]]$Nnode+Ntip(hp.cophylo$trees[[1]]),function(n,e) which(e==n),e=hp.cophylo$trees[[1]]$edge[,2]),frame="none",cex=2)

edgelabels.cophylo(hp.cophylo$trees[[2]]$node.label[2:hp.cophylo$trees[[2]]$Nnode],edge=sapply(2:hp.cophylo$trees[[2]]$Nnode+Ntip(hp.cophylo$trees[[2]]),function(n,e) which(e==n),e=hp.cophylo$trees[[2]]$edge[,2]),frame="none",which="right",cex=2)

dev.off()

pdf(file="tree1.pdf",width=50,height=80)

plotTree(tree1,color=NULL,fsize=0.5,ftype="reg",lwd=1,pts=FALSE)
nodelabels(text = tree1$node.label,frame = "n", cex=0.4)

dev.off()

pdf(file="tree2.pdf",width=50,height=80)

plotTree(tree2,color=NULL,fsize=0.5,ftype="reg",lwd=1,pts=FALSE)
nodelabels(text = tree1$node.label,frame = "n", cex=0.4)

dev.off()
