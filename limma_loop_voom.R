#!/usr/bin/env Rscript

#written for task 75.13

#usage:
#put the script somewhere in your PATH

#limma_loop_voom.R input.table metadata.table output_folder/ factors.list correction.factor 5

#'5'is the smallest allowable level size for each factor after missing value removal
#input data table is tab delimited with 'species' followed by sample names in first row, species name in first column and raw count data in each cell.
#factors list is a file of the factor column names from the metadata file to test as limma groups.. minus any correction factor variable

#missing data in metadata must be coded as an 'x', lowercase without the quotes.
#if you do not want a correction variable, then put empty quotes in the command ""

#N.B. Advice is to have any correction variable as a continuous variable, not a factor.

library(limma)
library(edgeR)
library(ggplot2)
library(data.table)
rm(list=ls()) # clear all

args = commandArgs(trailingOnly=TRUE)
setwd("./")
getwd()

otufile<-args[1]
metafile<-args[2]
outfolder<-args[3]
varslist<-args[4]
cor_var<- args[5]
min_count<-as.integer(args[6])
filt1<-args[7]

otu <- read.table(file = otufile, header=TRUE, sep="\t",row.names=1,check.names=FALSE,stringsAsFactors = F)

metadata <- read.table(file = metafile, header=TRUE, sep="\t", row.names=1,check.names=FALSE) 

varlist <- read.table(file = varslist,header=FALSE)

dgeoutfile <- DGEList(otu)#, group = testvar_)
dgeoutfileTMM <- edgeR::calcNormFactors(dgeoutfile, method = "TMM") 
norm_table <- cpm(dgeoutfileTMM)
write.table(norm_table,paste(outfolder,"tmm.tab",sep=""),sep="\t")
		

#remove samples with zero counts and corresponding metadata sample rows
sums<-colSums(otu)
zerosums<-sums!=0
otu<-otu[zerosums]
tmeta<-as.data.frame(t(metadata))
tmeta<-tmeta[zerosums]
metadata<-as.data.frame(t(tmeta))
metadata<-as.data.frame(apply(metadata,2,function(p)gsub('\\s+', '',p))) #remove whitepsace bug in transpose
#write.table(metadata,"meta.edit",sep="\t")
normdone<-0

n=0
for (i in varlist[,1]) {
  if (cor_var !=""){
	age_<-metadata[[cor_var]]}

  n=n+1
  
  j=as.character(varlist[n,])
  print (paste(n,j))

  testvar<-as.factor(metadata[[j]])
  meta<-metadata[[j]]
  
  #remove metadata NAs.. must be coded as 'x'#################################
  nas<-testvar!='x'
  
  otu2<-otu[,nas]
  testvar_<-testvar[nas]
  
  testvar_<-as.vector(testvar_)
  testvar_<-as.factor(testvar_)
  
  if (cor_var !=""){
      age_<-age_[nas]}
  
  ########################
  #check for all level size >min_count
  n2<-length(levels(testvar_))
  while (n2 >1) { #added test loop
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
  n2<-n2-1
  otu2<-otu2[,t2]
  if (cor_var !=""){
      age_<-age_[t2]}
  }
  
  ##ignore otus with no level with more than n non zeros
  otu3<-otu2[0,]
  removed<-0
  for (y in 1:nrow(otu2)){ 
  	otu4<-otu2[y,]
  	lvlgood<-0
  	for (z in levels(testvar_)){
 			lvl<-testvar_==z
 			sum1<-sum(otu4[lvl]!=0) #count  non zeros
  		if (sum1>=min_count){
  			lvlgood<-1} #the otu has at least one good level
  		}
  	
  	if (lvlgood==1){
  		otu3<-rbind(otu3,otu4)}
  	else{
  		removed<-removed+1}
  }
  print(paste(removed,"otus removed"))
  
  #repeat remove samples with library size zero
  sums1<-colSums(otu3)
	zerosums1<-sums1!=0
	otu3<-otu3[zerosums1]
	testvar_<-testvar_[zerosums1,drop=TRUE]
  if (cor_var!=""){
  	age_<-age_[zerosums1]
  }
  
  
  if (length(levels(testvar_)) <2 | cor_var==j |nrow(otu3)<2) {
    print (paste("only one level remaining in",j,"skipped"))
 
  } else {
  
    if (cor_var !=""){
   
      #age_<-as.factor(age_)
    	
      age_<-as.numeric(age_)
      
      design <- model.matrix(~0 + testvar_ + age_, data=metadata) #correction one factor
    }
    
    if (cor_var ==""){
      
    	#testvar can be continuous? No, because tests use contrasts only
    	
    	
      design <- model.matrix(~0 + testvar_ , data=metadata) #no correction
    }
    
    contdesign<-model.matrix(~0 + testvar_, data=metadata)
    
    dge <- DGEList(otu3, group = testvar_) 
    
	dgeTMM <- edgeR::calcNormFactors(dge, method = "TMM") 
	
    #keep <- filterByExpr(dge,group=testvar_)
    #dge<-dge[keep,]
    
    ## all combinations of groups
    combinations <- combn(rev(colnames(contdesign)), 2, function(x){paste(x, collapse = " - ")}) 
    cont1 <- makeContrasts(contrasts = combinations, levels=design)
    
    ## Running limma
    
    
    
	#dispersion plot - expect trend to be small
	#disp <- estimateDisp(dgeTMM, design = design, robust=TRUE)
    #disp[disp$prior.df==min(disp$prior.df),]
	#plotBCV(disp)
	
	v_OTU <- voom(dgeTMM, design = design, plot = FALSE)
    
	#if (n==1){
		#write.table(v_OTU$E,paste(outfolder,"voom.tab",sep=""),sep="\t")} #N.B. first var must have no missing values in order to get a complete table
	
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
    DT <- decideTests(fit)
    summary1<-summary(DT)
    
        cont_len=length(combinations)
    
    for (x1 in seq(1,cont_len)){ 
    
    print(x1)
	DT <- decideTests(fit)
    summary1<-summary(DT)
	
    top<-topTable(fit,n=Inf,coef=combinations[x1])
    
    
    varname<-gsub("testvar_","",combinations[x1])
    varname<-gsub(" ","",varname)
    
    
    write.table(top,paste(outfolder,j,"_",varname,'.limma',sep=""),sep="\t")
	
    }
    
	
  }
}
