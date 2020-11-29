#!/stornext/System/data/apps/R/R-3.5.2/lib64/R/bin/Rscript

#source("https://bioconductor.org/biocLite.R")
#biocLite("ALDEx2")

#devtools::install_github("Bioconductor-mirror/ALDEx2")

#install.packages("devtools") # if not already installed
#devtools::install_github("biomformat", "joey711")

library(biomformat)
library("ALDEx2")

args = commandArgs(trailingOnly=TRUE)

infile1<-"vir-met.biom"
outfile1<-"vir.aldex"
factor1<-"T1Dstatus"

dat <- read_biom(infile1)
otu_table <- as.data.frame(as.matrix(biom_data(dat)))
metadata <- sample_metadata(dat)
group<-metadata[factor1]
group<-group[,1]
tt <- aldex(otu_table, group,test = "t", mc.samples = 128)
aldex.plot(tt, type = "MW", cutoff = .05 / ncol(otu_table))
write.table(tt,outfile1,sep="\t")
