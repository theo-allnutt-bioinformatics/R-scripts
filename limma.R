#!/usr/bin/env Rscript

#install.packages("devtools") # if not already installed
#devtools::install_github("biomformat", "joey711")

library(biomformat)
library(phyloseq)
library(limma)
library(edgeR)
library(ggplot2)

args <- commandArgs(trailingOnly=TRUE)
setwd("./")

print (args[1])
print (args[2])

otu_table <- read.delim(file=args[1], header = TRUE, sep="\t", row.names=1)

metadata <- read.delim(file = args[2], header = TRUE, sep = "\t", quote="",row.names = 1) 

## A given design
metadata$parasite<-as.factor(metadata$parasite)
metadata$yob<-as.factor(metadata$yob)

design <- model.matrix(~ 0 + parasite + yob, data = metadata) # +DaysG + Nulliparous + BMI_conception

#design <- model.matrix(~T1Dstatus, data = metadata)
## Desired Contrast if fastor coded as 0/1 then add after name.. 'design will show what they are labeled as
cont <- makeContrasts(parasite = Parasite0 - Parasite1, levels = design)

## Running limma (example with T1Dstatus)
dge <- DGEList(otu_table, group = metadata$Parasite) #all samples must have at least one count

dgeTMM <- calcNormFactors(dge, method = "TMM")

v_OTU <- voom(dgeTMM, design = design, plot = TRUE)

norm_table <- cpm(dgeTMM)

write.table(norm_table,"tmm.tab",sep="\t")

## MAIN CHANGES INTRODUCED BY GORDON ##

## Voom pipeline with structural zeros --> fixing the residual
## df issue NOTE: a similar code would work for the
## limma-trend pipeline using the logCPM values

PoissonFit <- glmFit(dgeTMM, design, dispersion = 0, prior.count = 0)
StructuralZero <- (PoissonFit$fitted.values < 1e-08 & dgeTMM$counts < 1e-08)

v_OTU_NA <- v_OTU
v_OTU_NA$E[StructuralZero] <- NA

#corfit_NA <- duplicateCorrelation(v_OTU_NA, design, block = metadata$mothers)  # If you need to account for repeated measurements use a block colum = mother otherwise use SampleID

fit_NA <- lmFit(v_OTU_NA, design = design)#, block = metadata$mothers, correlation = corfit_NA$consensus.correlation)
fit <- lmFit(v_OTU, design = design)#, block = metadata$mothers, correlation = corfit_NA$consensus.correlation)

## Some of the elements of the fitted data are changed by the
## fitted data using NA instead of 0’s. According to Gordon,
## this is so not all the “0’s are being included in the
## regression”
fit$sigma <- fit_NA$sigma
fit$df.residual <- fit_NA$df.residual
fit$Amean <- fit_NA$Amean

fit <- contrasts.fit(fit, contrasts = cont)
fit <- eBayes(fit, robust = FALSE)
DT <- decideTests(fit)
summary(DT)

#write.fit(fit,DT,'fit.out')


top<-topTable(fit,n=Inf)

write.table(top,args[3],sep="\t")




