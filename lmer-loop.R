#!/usr/bin/env Rscript
#written for 75.19/

#reml with correction
#method from here https://ourcodingclub.github.io/tutorials/mixed-models/

library(ggplot2)
library(data.table)
library(lme4)

rm(list=ls()) # clear all

args = commandArgs(trailingOnly=TRUE)
setwd("./")
getwd()

otufile<-args[1] 
metafile<-args[2] 
outfolder<-args[3] 
varslist<-args[4]
corr1<- args[5]  
mincount<-args[6]

otu <- read.table(file = otufile, header=TRUE, sep="\t",row.names=1,check.names=FALSE,stringsAsFactors = F)

metadata <- read.table(file = metafile, header=TRUE, sep="\t", row.names=1,check.names=FALSE) 

varlist <- read.table(file = varslist,header=FALSE)


#remove samples with zero counts and corresponding metadata samples
#sums<-colSums(otu)
#zerosums<-sums!=0
#otu<-otu[zerosums]
#tmeta<-as.data.frame(t(metadata))
#tmeta<-tmeta[zerosums]
#metadata<-as.data.frame(tmeta) #t
#metadata<-as.data.frame(apply(metadata,2,function(p)gsub('\\s+', '',p))) #remove whitepsace bug in transpose
#write.table(metadata,"meta.edit",sep="\t")

#loop vars
n=0
for (i in varlist[,1]) {
  
	n=n+1
	j=as.character(varlist[n,])
	print (j)
  
	if (j!=corr1){
	
	otu_names<-c(rownames(otu)) 
	
	res1<-data.frame("deviance"= numeric(0), "Chisq"= numeric(0),"Pr(>Chisq)"=numeric(0))
	

	
	#loop otus
	#report test on each otu
	n1=0
	for (k in otu_names) {
		
		n1<-n1+1
		k<-otu_names[n1]
		#print(k)
	
			#remove data with missing values
		testvar_<-metadata[,j]
		data_<-otu
		nas<-testvar_!='x'   #remove missing values coded as 'x'
		data_<-data_[nas]
		testvar_<-testvar_[nas]
		testvar_<-as.vector(testvar_)
		if (corr1!=""){
			covar<-metadata[,corr1]
			covar<-covar[nas]}
		
		counts<-as.numeric(data_[k,])
		 ########################
  #check for all level size >min_count
	testvar2<-as.factor(testvar_)
  n2<-length(levels(testvar2))
  while (n2 >1) { #added test loop
  testcount<-min(table(testvar2))
  #print (testcount)
  #if (testcount<mincount){
    #print (paste("< ",mincount, " samples in some levels", testcount, "level excluded"))
  #}
  t2<-testvar2 
  levels(t2)[table(t2) < mincount]<-"TRUE" #########remove levels < min_count frequency
  t2 <- t2!="TRUE"
  testvar_<-testvar2[t2]
  testvar_<-as.vector(testvar_)
  testvar_<-as.factor(testvar_)
  n2<-n2-1
  counts<-counts[t2]
  if (corr1 !=""){
      covar<-covar[t2]}
  } #levels check while loop
  
  # loop level contrasts
  
  for (x1 in levels(testvar_)){
  
  
  if (length(levels(testvar_)) >1 ) {
    
  	######### end levels check run var if ok
		
		#check number of data points in model
		nonzero<-sum(counts>0, na.rm = TRUE)
		fvar<-as.factor(testvar_)
		nlev<-nlevels(fvar)
		
		#standardise variable
		testvar<-as.numeric(testvar_) 
		testvar <- scale(testvar, center = TRUE, scale = TRUE)
		if (corr1!=""){
				covar <- scale(covar, center = TRUE, scale = TRUE)}
	
		
		if (corr1=="" & nonzero >= mincount & nlev > 1){
		
			a.lm <- lm(counts ~ testvar) #, data = data1
			
			p<-as.data.frame(t(summary(a.lm)$coefficients[2,]),row.names=k)
			
			res1<-rbind(res1,p)
			}
		
		if (corr1!="" & nonzero >= mincount & nlev > 1){
			
			b.lm <- lm(counts ~ testvar + covar) #, data = data1
			
			p<-as.data.frame(t(summary(b.lm)$coefficients[2,]),row.names=k)
			
			res1<-rbind(res1,p)
			}
			
		else { 
			p<-as.data.frame(t(c(0,0,'NA','NA')),row.names=k)
			names(p)<-c("Estimate","Std. Error","t value","Pr(>|t|)")
			res1<-rbind(res1,p)
		}
		
  	}#levels check
		else{
			#print (paste("only one level remaining in",j,"skipped"))
			p<-as.data.frame(t(c(0,0,'NA','NA')),row.names=k)
			names(p)<-c("Estimate","Std. Error","t value","Pr(>|t|)")
			res1<-rbind(res1,p)
			}
		} #otu loop
		
		outfile<-paste(outfolder,j,".lmer",sep="")
		res1<-res1[order(res1$`Pr(>|t|)`),]
		
		#bh correction
		pa<-p.adjust(res1$`Pr(>|t|)`,method="fdr")
		
		res1<-cbind(res1,"fdr"=pa)
		
		write.table(res1,outfile,append=FALSE,quote=FALSE,sep="\t")
		
	} #check testvar<>covar
	
	} #var loop

