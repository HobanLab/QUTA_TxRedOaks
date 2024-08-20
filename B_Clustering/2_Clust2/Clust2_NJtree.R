# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% QUTA & TEXAS RED OAKS : BUILDING NJ TREE %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# This is an exploratory script, for building a NJ trees of QUTA and Texas Red Oak samples
# dataset. Different datasets are grouped into different sections

library(adegenet)
library(ape)
library(Rphylip)
library(seqinr)
library(rdiversity)
library(vcfR)

# Folder containing different outputs
QUTA_TRO_wd <- 
  '/RAID1/QUTA_TX_RedOaks/Genotyping/SNP_Calling/Output/'
setwd(QUTA_TRO_wd)

# SE-READS ONLY, FIRST CLUSTERING RUN ----
seOnly_Clust1_folder <- paste0(QUTA_TRO_wd,'All_seReadsOnly/pop_R97/')
# Read in genind file
QUTA_TRO_gen <- 
  read.genepop(paste0(seOnly_Clust1_folder, "populations.snps.gen"))
# Convert genind file to data.frame
QUTA_TRO_df <- genind2df(QUTA_TRO_gen)
# Convert data.frame to genlight object
QUTA_TRO_genLight <- as.genlight(QUTA_TRO_df)
# Build a genetic distance matrix from the genlight object
QUTA_TRO_distMat <- dist(QUTA_TRO_genLight)
# Build a NJ tree from the genetic distance matrix, and plot
QUTA_TRO_njTree <- nj(QUTA_TRO_distMat)
plot(QUTA_TRO_njTree, "phylo")

# SE-READS ONLY, SECOND CLUSTERING RUN, WHITELISTED LOCI ----
seOnly_Clust2_folder <- paste0(QUTA_TRO_wd,'Clust2_seOnly/pop_R98_WL10k/')
# Read in genind file
QUTA_TRO_gen <- 
  read.genepop(paste0(seOnly_Clust2_folder, "populations.snps.gen"))
# Convert genind file to data.frame
QUTA_TRO_df <- genind2df(QUTA_TRO_gen)
# Convert data.frame to genlight object
QUTA_TRO_genLight <- as.genlight(QUTA_TRO_df)
# Build a genetic distance matrix from the genlight object
QUTA_TRO_distMat <- dist(QUTA_TRO_genLight)
# Build a NJ tree from the genetic distance matrix
QUTA_TRO_njTree <- nj(QUTA_TRO_distMat)
# Rename the phylognetic tree, to make more discernible
QUTA_TRO_njTree$tip.label <- sub('SYST-MOR-000','',QUTA_TRO_njTree$tip.label)
QUTA_TRO_njTree$tip.label <- sub('OAK-MOR-00','',QUTA_TRO_njTree$tip.label)
QUTA_TRO_njTree$tip.label <- sub('OAK-MOR-000','',QUTA_TRO_njTree$tip.label)
QUTA_TRO_njTree$tip.label <- sub('OAK-UMN-000','',QUTA_TRO_njTree$tip.label)
# Plot tree
plot(QUTA_TRO_njTree, "phylo")

# SE-READS ONLY, SECOND CLUSTERING RUN, ALL LOCI ----
seOnly_Clust2_folder <- paste0(QUTA_TRO_wd,'Clust2_seOnly/pop_R98/')
# Read in genind file
QUTA_TRO_gen <- 
  read.genepop(paste0(seOnly_Clust2_folder, "populations.snps.gen"))
# Convert genind file to data.frame
QUTA_TRO_df <- genind2df(QUTA_TRO_gen)
# Convert data.frame to genlight object
QUTA_TRO_genLight <- as.genlight(QUTA_TRO_df)
# Build a genetic distance matrix from the genlight object
QUTA_TRO_distMat <- dist(QUTA_TRO_genLight)
# Build a NJ tree from the genetic distance matrix
QUTA_TRO_njTree <- nj(QUTA_TRO_distMat)
# Rename the phylognetic tree, to make more discernible
QUTA_TRO_njTree$tip.label <- sub('SYST-MOR-000','',QUTA_TRO_njTree$tip.label)
QUTA_TRO_njTree$tip.label <- sub('OAK-MOR-00','',QUTA_TRO_njTree$tip.label)
QUTA_TRO_njTree$tip.label <- sub('OAK-MOR-000','',QUTA_TRO_njTree$tip.label)
QUTA_TRO_njTree$tip.label <- sub('OAK-UMN-000','',QUTA_TRO_njTree$tip.label)
# Plot tree
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
