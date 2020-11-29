#!/stornext/System/data/apps/R/R-3.5.2/lib64/R/bin/Rscript

#install.packages("devtools") # if not already installed
#devtools::install_github("biomformat", "joey711")

args = commandArgs(trailingOnly=TRUE)

library(biomformat)
library(phyloseq)
library("limma", lib.loc="/stornext/System/data/apps/R/R-3.5.2/lib64/R/library")
library("edgeR", lib.loc="/stornext/System/data/apps/R/R-3.5.2/lib64/R/library")


otu_table_met <- import_biom("virus_refseq_s.biom")     #args[1])


#otu_table <- merge_samples(Otus_table_C, "SampleID", fun = sum)
#Otus_table_MC <- prune_taxa(taxa_sums(Otus_table_MC) >= 1, Otus_table_MC) #remove otus with fewer than given number accross all samples
metadata <- read.delim(file = "mapping2.txt", header = TRUE, sep = "\t", quote="",row.names = 1) # '#'must be removed from samplid

Otus_table_MC <- phyloseq(otu_table(otu_table_met), tax_table(otu_table_met), sample_data(metadata))

Counts <- as.data.frame(otu_table(Otus_table_MC))  # Using the OTU table that wasn't normalised. You can use the rarefied one in case you had the need to normalise
# Metadata
Meta_df <- as.matrix(sample_data(Otus_table_MC))
Meta_df <- as.data.frame(Meta_df,stringsAsFactors = FALSE)
## A given design
design <- model.matrix(~0 + T1Dstatus + DaysG + Nulliparous + Age_LMP + BMI_conception + HLA.6DRML, data = Meta_df)

## Desired Contrast, if you need it
cont <- makeContrasts(nonT1D_vs_T1D = T1Dstatus0 - T1Dstatus1, levels = design)

#working to here########################################


## Running limma (example with T1Dstatus)
dge <- DGEList(Counts, group = Meta_df$T1Dstatus)
dgeTMM <- calcNormFactors(dge, method = "TMM")
v_OTU <- voom(dgeTMM, design = design, plot = FALSE)

## MAIN CHANGES INTRODUCED BY GORDON ##

## Voom pipeline with structural zeros --> fixing the residual
## df issue NOTE: a similar code would work for the
## limma-trend pipeline using the logCPM values

PoissonFit <- glmFit(dgeTMM, design, dispersion = 0, prior.count = 0)
StructuralZero <- (PoissonFit$fitted.values < 1e-08 & dgeTMM$counts < 
                     1e-08)

v_OTU_NA <- v_OTU
v_OTU_NA$E[StructuralZero] <- NA

corfit_NA <- duplicateCorrelation(v_OTU_NA, design, block = Meta_df$motherid)  # If you need to account for repeated measurements

fit_NA <- lmFit(v_OTU_NA, design = design, block = Meta_df$motherid, correlation = corfit_NA$consensus.correlation)
fit <- lmFit(v_OTU, design = design, block = Meta_df$motherid, correlation = corfit_NA$consensus.correlation)

## Some of the elements of the fitted data are changed by the
## fitted data using NA instead of 0’s. According to Gordon,
## this is so not all the “0’s are being included in the
## regression”
fit$sigma <- fit_NA$sigma
fit$df.residual <- fit_NA$df.residual
fit$Amean <- fit_NA$Amean

fit <- contrasts.fit(fit, contrasts = cont)
fit <- eBayes(fit, robust = FALSE)
DT <- decideTests(fit)
summary(DT)










