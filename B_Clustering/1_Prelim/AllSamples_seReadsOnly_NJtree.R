# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% QUTA & TEXAS RED OAKS : BUILDING NJ TREE %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# This is an exploratory script, for building a NJ tree of QUTA and Texas Red Oak samples
# from single-read data only, using a missing data filer of R97

library(adegenet)
library(ape)
library(Rphylip)
library(seqinr)
library(rdiversity)
library(vcfR)

QUTA_TRO_wd <- '/RAID1/QUTA_TX_RedOaks/Genotyping/SNP_Calling/Output/All_seReadsOnly/pop_R97/'
setwd(QUTA_TRO_wd)

# USING ADEGENET ----
# Read in genind file
QUTA_TRO_gen <- read.genepop("populations.snps.gen")
# Convert genind file to data.frame
QUTA_TRO_df <- genind2df(QUTA_TRO_gen)
# Convert data.frame to genlight object
QUTA_TRO_genLight <- as.genlight(QUTA_TRO_df)
# Build a genetic distance matrix from the genlight object
QUTA_TRO_distMat <- dist(QUTA_TRO_genLight)
# Build a NJ tree from the genetic distance matrix, and plot
QUTA_TRO_njTree <- nj(QUTA_TRO_distMat)
plot(QUTA_TRO_njTree, "phylo")

# OTHER APPROACHES (UNSUCCESSFUL) ----
# Read in the phylip file (as a data.frame)...
QUTA_TRO_njTree <- read.phylip("populations.all.phylip")

# seqinr
QUTA_TRO_njTree <- read.alignment("populations.all.phylip", format = "phylip")
# Error: Error in read.alignment("populations.all.phylip", format = "phylip"): 
# File populations.all.phylip is not readable

# rdiversity
QUTA_TRO_vcf <- read.vcfR("populations.snps.vcf")
QUTA_TRO_distMat <- gen2dist(QUTA_TRO_vcf, biallelic = TRUE)
# Error: Error in h(simpleError(msg, call)): 
# error in evaluating the argument 'i' in selecting a method for function '[': argument of length 0
