#!/usr/bin/env Rscript

library(biomformat)
library(phyloseq)
library(limma)
library(edgeR)
library(ggplot2)

rm(list=ls()) # clear all
args <- commandArgs(trailingOnly=TRUE)

## OTU table
OTU_Met <- read.table(file=args[1], header=TRUE, sep="\t", row.names=1,check.names=FALSE)
## Metadata



## DIFFERENTIAL ABUNDANCE with Limma
#Come back to data frame
Counts2 <- as.data.frame(otu_table(OTU_Met_F))


### 1) Limma with design 1
dge <- DGEList(Counts2)
dgeTMM <- edgeR::calcNormFactors(dge, method = "TMM")

norm_table <- cpm(dgeTMM)

write.table(norm_table,paste(args[1],".tmm"),sep="\t")

