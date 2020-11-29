#!/usr/bin/env Rscript
rm(list=ls()) # clear all
args = commandArgs(trailingOnly=TRUE)
setwd("./")
getwd()

# This runs SourceTracker on the original "contamination" data set
# (data included in 'data' folder)



# load sample metadata

otufile<-args[1]
mappingfile<-args[2]
outfile<-args[3]

metadata <- read.table(mappingfile,sep='\t',h=T,row.names=1,check=F,comment='')

# load OTU table
# This 'read.table' command is designed for a 
# QIIME-formatted OTU table.
# namely, the first line begins with a '#' sign
# and actually _is_ a comment; the second line
# begins with a '#' sign but is actually the header
otus <- read.table(otufile,sep='\t', header=T,row.names=1,check=F,skip=1,comment='')
otus <- t(as.matrix(otus))

# extract only those samples in common between the two tables
common.sample.ids <- intersect(rownames(metadata), rownames(otus))
otus <- otus[common.sample.ids,]
metadata <- metadata[common.sample.ids,]
# double-check that the mapping file and otu table
# had overlapping samples
if(length(common.sample.ids) <= 1) {
    message <- paste(sprintf('Error: there are %d sample ids in common '),
                    'between the metadata file and data table')
    stop(message)
}

# extract the source environments and source/sink indices
train.ix <- which(metadata$SourceSink=='source')
test.ix <- which(metadata$SourceSink=='sink')
envs <- metadata$Env
if(is.element('Description',colnames(metadata))) desc <- metadata$Description


# load SourceTracker package
source('/stornext/HPCScratch/home/allnutt.t/bin/sourcetracker-1.0.1/src/SourceTracker.r')

# tune the alpha values using cross-validation (this is slow!)
# tune.results <- tune.st(otus[train.ix,], envs[train.ix])
# alpha1 <- tune.results$best.alpha1
# alpha2 <- tune.results$best.alpha2
# note: to skip tuning, run this instead:
alpha1 <- alpha2 <- 0.001

# train SourceTracker object on training data
st <- sourcetracker(otus[train.ix,], envs[train.ix])

# Estimate source proportions in test data
results <- predict(st,otus[test.ix,], alpha1=alpha1, alpha2=alpha2)

# Estimate leave-one-out source proportions in training data 
results.train <- predict(st, alpha1=alpha1, alpha2=alpha2)

# plot results
 

labels <- sprintf('%s %s', envs,desc)
#pie1<-plot(results, labels[test.ix], type='pie')

# other plotting functions
#bar1<-plot(results, labels[test.ix], type='bar')
#dist1<-plot(results, labels[test.ix], type='dist')
#pietrain<-plot(results.train, labels[train.ix], type='pie')
#bartrain<-plot(results.train, labels[train.ix], type='bar')
#disttrain<-plot(results.train, labels[train.ix], type='dist')

# plot results with legend
#legendplot<-plot(results, labels[test.ix], type='pie', include.legend=TRUE, env.colors=c('#47697E','#5B7444','#CC6666','#79BEDB','#885588'))

#pdf(paste(outfile,"_","pie",".pdf",sep=""))
			#print(pie1)     
#pdf(paste(outfile,"_","bar",".pdf",sep=""))
			#print(bar1) 
#pdf(paste(outfile,"_","dist",".pdf",sep=""))
			#print(dist1)
#pdf(paste(outfile,"_","pie-train",".pdf",sep=""))
			#print(pietrain)
#pdf(paste(outfile,"_","bar-train",".pdf",sep=""))
			#print(bartrain)
#pdf(paste(outfile,"_","dist-train",".pdf",sep=""))
			#print(disttrain)
#pdf(paste(outfile,"_","legend",".pdf",sep=""))
			#print(legendplot)

#dev.off()

write.table(results$proportions,paste(outfile,"prop.tab",sep=""),sep="\t")
write.table(results$proportions_sd,paste(outfile,"prop_sd.tab",sep=""),sep="\t")

write.table(results.train$proportions,paste(outfile,"train-prop.tab",sep=""),sep="\t")
write.table(results.train$proportions_sd,paste(outfile,"train-prop_sd.tab",sep=""),sep="\t")


