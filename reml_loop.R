#!/usr/bin/env Rscript
#written for 75.1/alpha_tests

#reml with correction
#method from here https://ourcodingclub.github.io/tutorials/mixed-models/

library(ggplot2)
library(data.table)
library(lme4)

rm(list=ls()) # clear all

args = commandArgs(trailingOnly=TRUE)
setwd("./")
getwd()

otufile<-args[1] #"/stornext/HPCScratch/home/allnutt.t/d/075_ELCHO_metagenome/75.1metaphlan2/alpha_tests/alphadiv/s.tab" 
metafile<-args[2] #"/stornext/HPCScratch/home/allnutt.t/d/075_ELCHO_metagenome/75.1metaphlan2/alpha_tests/mapping_reml.txt"
outfile<-args[3] #"/stornext/HPCScratch/home/allnutt.t/d/075_ELCHO_metagenome/75.1metaphlan2/alpha_tests/reml.txt"
varlist<-args[4]
corr1<-args[5] #"age_qtrs"
 #"/stornext/HPCScratch/home/allnutt.t/d/075_ELCHO_metagenome/75.1metaphlan2/alpha_tests/vars.txt" 

otu <- read.table(file = otufile, header=TRUE, sep="\t",row.names=1,check.names=FALSE,stringsAsFactors = F)

metadata <- read.table(file = metafile, header=TRUE, sep="\t", row.names=1,check.names=FALSE) 

varlist <- read.table(file = varlist,header=FALSE)

write.table(paste("fixed","richness","shannon2","simpson",sep="\t"),outfile,col.names=FALSE,row.names=FALSE,quote=FALSE)

#loop vars
n=0
for (i in varlist[,1]) {
  
	n=n+1
	j=as.character(varlist[n,])
	print (j)
  
	alphadivs<-c(colnames(otu))
	#loop alpha measures
	res1<-list()
	for (k in alphadivs) {
		
		print(k)
		data1<-cbind(otu,metadata)

		testscore<-data1[[k]]
		testvar<-data1[[j]] #age_months #days_on_AB
		
		if (corr1!=""){
		covar<-data1[[corr1]]}

		#not sure if we need to scale variable
		#testvar <- scale(testvar, center = TRUE, scale = TRUE)

		#remove NAs.. must be coded as 'x'########
		nas<-testvar!='x'
		testscore<-testscore[nas]
		testvar<-testvar[nas]
		if (corr1!=""){
		covar<-covar[nas]}
		
		data1<-data1[nas,]
		
		if (corr1!=""){
		a.lm <- lmer(testscore ~ testvar +(1|covar), data = data1,REML=FALSE)
		b.lm <- lmer(testscore ~ 1 +(1|covar), data = data1,REML=FALSE)
		ab<-anova(a.lm,b.lm)
		p<-ab$`Pr(>Chisq)`[2]
		res1[[k]]<-p
		} #covar if1
		
		if (corr1==""){
		a.lm <- lm(testscore ~ testvar, data = data1)
		p<-summary(a.lm)$coefficients[2,4]
		
		res1[[k]]<-p
		}#covar if2
		
		if (p<=0.05){
			yvar<-testscore
			xvar<-testvar
			nas<-xvar!="x"
			xvar<-xvar[nas]
			yvar<-yvar[nas]
			xvar<-as.factor(xvar)
			data1<-data1[nas,]
			p<-round(p,digits=5)
			
			aplot<-ggplot(data1, aes(x=xvar, y=yvar)) + geom_boxplot() + geom_jitter(width=0.1) +theme(axis.text.x = element_text(angle = 90),panel.background = element_rect(fill = "white", colour = "black",size = 1, linetype = "solid"),panel.grid.major = element_line(size = 0.5, linetype = 'solid',colour = "grey"),panel.grid.minor = element_line(size = 0.25, linetype = 'solid',colour = "grey"))+labs(x=j,y=paste(k,"p =",p))
			
			pdf(paste(outfile,"_",j,".pdf",sep=""))
			print(aplot)     
			dev.off() 

			}#plot if
		
		
		
		
		
		} #alpha div loop
	output1<-paste(j,res1[["richness"]],res1[["shannon_2"]],res1[["simpson"]],sep="\t")
	
	write.table(output1,outfile,append=TRUE,col.names=FALSE,row.names=FALSE,quote=FALSE)
	} #var loop


#plot(a.lm) 
#qqnorm(resid(a.lm))
#qqline(resid(a.lm))





