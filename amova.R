#!/usr/bin/Rscript

library(pegas)

args = commandArgs(trailingOnly=TRUE)
reps=as.numeric(args[3])
all.data <- read.table(file=args[1],sep="\t",header=TRUE)

p<-all.data$factor

all.data <- all.data[,-1]
d=as.matrix(all.data)

am <- vector("list", 15)
amova1<- amova(d~p,nperm=reps,is.squared=TRUE)
counter <- 1
n=as.numeric(args[4])
for(i in 1:(n-1)){
   for(j in (i+1):n){
#if(i != j){
p1 = (p == unique(p)[i])
p2 = (p == unique(p)[j])
a=d[p1 | p2, p1 | p2]
b=p[p1 | p2]
print (a)
print (c(counter,b))
am[[counter]] <- amova(a ~ b,nperm=reps, is.squared = TRUE)
print (am[counter])
counter <- counter + 1
}}
pvalues=unlist(lapply(am, function(i) i$varcomp$P.value[1]))
phist=unlist(lapply(am, function(i) i$varcomp$sigma2[1]/(i$varcomp$sigma2[1]+i$varcomp$sigma2[2])))
print (phist)
print (pvalues)
#force into table

phimat<-matrix(,nrow=n,ncol=n)
pmat<-matrix(,nrow=n,ncol=n)
cnt=1
for (x in 1:(n)){
for (y in x:(n)){
print (c(cnt,x,y))
if (x==y){
phimat[x,y]=0
pmat[x,y]=0
}else{
phimat[x,y]=phist[cnt]
pmat[x,y]=pvalues[cnt]
cnt=cnt+1
}
}}
print (amova1)
print (phimat)
print (pmat)
phimat<-t(phimat)
pmat<-t(pmat)
a1="AMOVA RESULT\n"
a2="PAIRWISE PHIST\n"
a3="\nPAIRWISE P-VALUES\n"

write(a1,file=args[2],ncolumns=1,append = FALSE,sep="\t")
capture.output(amova1, file = args[2],append=TRUE)
write(a2,file=args[2],ncolumns=1,append = TRUE,sep="\t")
write(levels(p),file=args[2],ncolumns=n,append = TRUE,sep="\t")
write(phimat,file=args[2],ncolumns=n,append = TRUE,sep="\t")
write(a3,file=args[2],ncolumns=1,append = TRUE,sep="\t")
write(levels(p),file=args[2],ncolumns=n,append = TRUE,sep="\t")
write(pmat,file=args[2],ncolumns=n,append = TRUE,sep="\t")
