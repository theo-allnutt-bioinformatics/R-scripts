#!/usr/bin/env Rscript

#written for task 75.13

#usage:
#put the script somewhere in your PATH
#limma_loop_voom.R input.table metadata.table output_folder/ factors.list correction.factor 5
#'5'is the smallest allowable level size for each factor
#input data table is tab delimited with 'species' followed by sample names in first row, species name in first column and raw count data in each cell.
#factors list is a file of the factor column names from the metadata to test as limma groups.. minus any correction factor variable


#install.packages("devtools") # if not already installed
#devtools::install_github("biomformat", "joey711")

#library(biomformat)

library(limma)
library(edgeR)
library(ggplot2)
library(data.table)
rm(list=ls()) # clear all

args = commandArgs(trailingOnly=TRUE)
setwd("/stornext/HPCScratch/home/allnutt.t/d/075_ELCHO_metagenome/75.13_unified_analysis/burst/")
getwd()

otufile<-"/stornext/HPCScratch/home/allnutt.t/d/075_ELCHO_metagenome/75.13_unified_analysis/burst/euk.tab"
metafile<-"/stornext/HPCScratch/home/allnutt.t/d/075_ELCHO_metagenome/75.13_unified_analysis/burst/mapping-euk.txt"
outfolder<-"/stornext/HPCScratch/home/allnutt.t/d/075_ELCHO_metagenome/75.13_unified_analysis/burst/limma-test/"
cor_var<- "age_months"
min_count<-as.integer(5)

otu <- read.table(file = otufile, header=TRUE, sep="\t",row.names=1,check.names=FALSE,stringsAsFactors = F)

metadata <- read.table(file = metafile, header=TRUE, sep="\t", row.names=1,check.names=FALSE) 

varlist <- read.table(file = "vars.txt",header=FALSE)


n=0
for (i in varlist[,1]) {
  
  n=n+1
  
  j=as.character(varlist[n,])
  print (j)
  
 
  testvar<-as.factor(metadata[[j]])
  meta<-metadata[[j]]
  
  #remove NAs.. must be coded as 'x'#################################################
  
  nas<-testvar!='x'
  
  otu2<-otu[,nas]
  testvar_<-testvar[nas]
  #age_<-age[nas]
  
  testvar_<-as.vector(testvar_)
  testvar_<-as.factor(testvar_)
  #age_<-as.factor(age_)
  
  #test for minim
  
  ########################
  #check for minimum level size >min_count
  testcount<-min(table(testvar_))
  #print (testcount)
  if (testcount<min_count){
    print (paste("< ",min_count, " samples in some levels", testcount, "level excluded"))
  }
  t2<-testvar_ 
  levels(t2)[table(t2) < min_count]<-"TRUE" #########remove levels < min_count frequency
  t2 <- t2!="TRUE"
  testvar_<-testvar_[t2]
  testvar_<-as.vector(testvar_)
  testvar_<-as.factor(testvar_)
  
  
  otu2<-otu2[,t2]
  
  
  if (length(levels(testvar_)) <2) {
    print (paste("only one level remaining in",j,"skipped"))
 
  } else {
    
    if (cor_var !=""){
      
      age_<-metadata[[cor_var]]
      
      age_<-age_[nas]
      age_<-age_[t2]
      #age_<-as.factor(age_)
      
      design <- model.matrix(~0 + testvar_ + age_, data=metadata) #correction one factor
    }
    
    if (cor_var ==""){
      
      design <- model.matrix(~0 + testvar_ , data=metadata) #no correction
    }
    
    contdesign<-model.matrix(~0 + testvar_, data=metadata)
    
    dge <- DGEList(otu2, group = testvar_) 
    
    ## all combinations of groups
    combinations <- combn(rev(colnames(contdesign)), 2, function(x){paste(x, collapse = " - ")}) 
    cont1 <- makeContrasts(contrasts = combinations, levels=design)
    
    ## Running limma
    
    dgeTMM <- edgeR::calcNormFactors(dge, method = "TMM") 
     if (n==1){
		norm_table <- cpm(dgeTMM)

		write.table(norm_table,paste(outfolder,"limma.tmm",sep=""),sep="\t")} #N.B. first var must have no missing values in order to get a complete table
    
	#dispersion plot - expect trend to be small
	#disp <- estimateDisp(dgeTMM, design = design, robust=TRUE)
    #disp[disp$prior.df==min(disp$prior.df),]
	#plotBCV(disp)
	
	v_OTU <- voom(dgeTMM, design = design, plot = FALSE)
    
	if (n==1){
		write.table(v_OTU$E,paste(outfolder,"limma.voom",sep=""),sep="\t")} #N.B. first var must have no missing values in order to get a complete table
	
    ## MAIN CHANGES INTRODUCED BY GORDON SMITH##
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
    #DT <- decideTests(fit)
    #summary1<-summary(DT)
    
    cont_len=length(combinations)
    
    for (x1 in seq(1,cont_len)){ 
    
    print(combinations[x1])
    	
    cont1 <- makeContrasts(contrasts = combinations[x1], levels=design)
    	
    DT <- decideTests(fit)
    summary1<-summary(DT)
    
    top<-topTable(fit,n=Inf,coef=combinations[x1])
    
    
    write.table(j,paste(outfolder,'limma-all-summary.txt',sep=""),sep="\t",append=TRUE)
    write.table(summary1,paste(outfolder,'limma-all-summary.txt',sep=""),sep="\t",append=TRUE)
    
    write.table(top,paste(outfolder,j,'.limma',sep=""),sep="\t")
    }
  }
}

