#!/usr/bin/Rscript


library("RnBeads")

args = commandArgs(trailingOnly=TRUE)

#options(fftempdir="/scratch/mandrea_data/temp")

rnb.options(disk.dump.big.matrices=FALSE)

logger.start(fname=NA)
parallel.isEnabled()
num.cores <- 8
parallel.setup(num.cores)
parallel.isEnabled()
data.dir <- args[1]
idat.dir <- file.path(data.dir, "idat")
sample.annotation <- file.path(data.dir, args[2])
analysis.dir <- args[3]
report.dir <- file.path(analysis.dir, "reports")
rnb.options(filtering.sex.chromosomes.removal=TRUE, identifiers.column="Sample_ID")
rnb.run.analysis(dir.reports=report.dir, sample.sheet=sample.annotation, data.dir=idat.dir, data.type="infinium.idat.dir")
