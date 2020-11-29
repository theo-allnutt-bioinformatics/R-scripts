#!/usr/bin/env Rscript
#written for 75.1/alpha_tests

#reml with correction

library(ggplot2)
library(data.table)

rm(list=ls()) # clear all

args = commandArgs(trailingOnly=TRUE)
setwd("/stornext/HPCScratch/home/allnutt.t/d/075_ELCHO_metagenome/75.1metaphlan2/alpha_tests/")
getwd()

otufile<-"/stornext/HPCScratch/home/allnutt.t/d/075_ELCHO_metagenome/75.1metaphlan2/alpha_tests/alphadiv/s.tab"
metafile<-"/stornext/HPCScratch/home/allnutt.t/d/075_ELCHO_metagenome/75.1metaphlan2/alpha_tests/mapping_reml.txt"
outfolder<-"/stornext/HPCScratch/home/allnutt.t/d/075_ELCHO_metagenome/75.1metaphlan2/alpha_tests/alphadiv/"
covar<-"age_months"
min_count<-as.integer(5)
varfile<-"vars.txt"
otu <- read.table(file = otufile, header=TRUE, sep="\t",row.names=1,check.names=FALSE,stringsAsFactors = F)

metadata <- read.table(file = metafile, header=TRUE, sep="\t", row.names=1,check.names=FALSE) 

#loop vars

#loop alpha measures

data1<-cbind(otu,metadata)

testscore<-data1$richness
testvar<-data1$Age_complementary_food_started #age_months #days_on_AB
covar<-data1[[covar]]

#not sure if we need to scale variable
#testvar <- scale(testvar, center = TRUE, scale = TRUE)

#remove NAs.. must be coded as 'x'#################################################
  nas<-testvar!='x'
  testscore<-testscore[nas]
  testvar<-testvar[nas]
  covar<-covar[nas]
  data1<-data1[nas,]

alpha.lm <- lm(testscore ~ testvar +covar, data = data1)
summary(alpha.lm)

#(prelim_plot <- ggplot(data1, aes(x = testvar, y = testscore)) + geom_point() + geom_smooth(method = "lm"))
#residuals - should be flat
#plot(basic.lm, which = 1) 
#plot(basic.lm, which = 2)




