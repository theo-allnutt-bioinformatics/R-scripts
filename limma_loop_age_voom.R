#!/usr/bin/env Rscript

#install.packages("devtools") # if not already installed
#devtools::install_github("biomformat", "joey711")

#library(biomformat)
library(phyloseq)
library(limma)
library(edgeR)
library(ggplot2)

rm(list=ls()) # clear all

args = commandArgs(trailingOnly=TRUE)
setwd("./")
getwd()

otufile<-args[1]
metafile<-args[2]
outfolder<-args[3]

otu <- read.table(file = args[1], header=TRUE, sep="\t",row.names=1,check.names=FALSE,stringsAsFactors = F)

metadata <- read.table(file = args[2], header=TRUE, sep="\t", row.names=1,check.names=FALSE) 

varlist <- read.table(file = args[4],header=FALSE)

n=0
for (i in varlist[,1]) {

n=n+1

j=as.character(varlist[n,])
print (j)

testvar<-as.factor(metadata[[j]])
meta<-metadata[[j]]

age<-metadata$age_yrs

#remove NAs.. must be coded as 'x'

nas<-testvar!='x'

otu2<-otu[,nas]
testvar_<-testvar[nas]
age_<-age[nas]

testvar_<-as.vector(testvar_)
testvar_<-as.factor(testvar_)
age_<-as.factor(age_)

########################
#check for minimum level size >3
testcount<-min(table(testvar_))
if (testcount<3){
  print (j,"too few samples, min=",testcount)
  print (table(testvar_))
  break
}

design <- model.matrix(~0 + testvar_ + age_, data=metadata)
contdesign<-model.matrix(~0 + testvar_, data=metadata)

dge <- DGEList(otu2, group = testvar_) 

## all combinations of groups
combinations <- combn(rev(colnames(contdesign)), 2, function(x){paste(x, collapse = " - ")}) 
cont1 <- makeContrasts(contrasts = combinations, levels=design)

## Running limma

dgeTMM <- edgeR::calcNormFactors(dge, method = "TMM") 

v_OTU <- voom(dgeTMM, design = design, plot = TRUE)

## MAIN CHANGES INTRODUCED BY GORDON ##
## Voom pipeline with structural zeros --> fixing the residual
## df issue NOTE: a similar code would work for the
## limma-trend pipeline using the logCPM values

PoissonFit <- glmFit(dgeTMM, design, dispersion = 0, prior.count = 0)
StructuralZero <- (PoissonFit$fitted.values < 1e-08 & dgeTMM$counts < 1e-08)

v_OTU_NA <- v_OTU
v_OTU_NA$E[StructuralZero] <- NA

#corfit_NA <- duplicateCorrelation(v_OTU_NA, design, block = metadata$motherid)# If you need to account for repeated measurements use a block colum = mother otherwise use SampleID

fit_NA <- lmFit(v_OTU_NA, design = design)#, block = metadata$motherid, correlation = corfit_NA$consensus.correlation)
fit <- lmFit(v_OTU, design = design)#, block = metadata$motherid, correlation = corfit_NA$consensus.correlation)

## Some of the elements of the fitted data are changed by the
## fitted data using NA instead of 0’s. According to Gordon,
## this is so not all the “0’s are being included in the
## regression”
fit$sigma <- fit_NA$sigma
fit$df.residual <- fit_NA$df.residual
fit$Amean <- fit_NA$Amean

fit <- contrasts.fit(fit, contrasts = cont1)
fit <- eBayes(fit, robust = FALSE)
DT <- decideTests(fit)
summary1<-summary(DT)

top<-topTable(fit,n=Inf,coef=combinations)

#write.table(summary1,paste(j,'-summary.limma'),sep="\t")
write.table(j,paste(outfolder,'limma-all-summary.txt'),sep="\t",append=TRUE)
write.table(summary1,paste(outfolder,'limma-all-summary.txt'),sep="\t",append=TRUE)

write.table(top,paste(outfolder,j,'.limma'),sep="\t")

}


