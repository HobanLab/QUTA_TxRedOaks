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

QUTA_TRO_wd <- '/RAID1/QUTA_TX_RedOaks/Genotyping/SNP_Calling/Output/Clust3_seOnly/pop_R98/'
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
# Build a NJ tree from the genetic distance matrix
QUTA_TRO_njTree <- nj(QUTA_TRO_distMat)
# Rename the phylognetic tree, to make more discernible
QUTA_TRO_njTree$tip.label <- sub('SYST-MOR-000','',QUTA_TRO_njTree$tip.label)
QUTA_TRO_njTree$tip.label <- sub('OAK-MOR-00','',QUTA_TRO_njTree$tip.label)
QUTA_TRO_njTree$tip.label <- sub('OAK-MOR-000','',QUTA_TRO_njTree$tip.label)
QUTA_TRO_njTree$tip.label <- sub('OAK-UMN-000','',QUTA_TRO_njTree$tip.label)
# Plot
par(mar=c(1,1,2,1)+0.1)
plot(QUTA_TRO_njTree, 'phylo', main='Clust3 NJ Tree: Abbreviated Sample Names')

# Read in a CSV listing sample names in the order indicated in QUTA_TRO_njTree$tip.label
clustNames <- read.csv2('/home/akoontz/Documents/QUTA_TxRedOaks/Code/clustIDs_withQUTA_numOrdered.csv', 
                        header = TRUE, sep = ",")[,2]
# Rename and plot
QUTA_TRO_njTree$tip.label <- clustNames
par(mar=c(1,1,2,1)+0.1)
plot(QUTA_TRO_njTree, 'phylo', main='Clust3 NJ Tree: K4 Cluster Identities')
