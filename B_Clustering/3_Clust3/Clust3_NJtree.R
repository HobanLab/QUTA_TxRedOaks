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

# CLUST3 %%%% ----
QUTA_Clust3_wd <- '/RAID1/QUTA_TX_RedOaks/Genotyping/SNP_Calling/Output/Clust3_seOnly/pop_R98/'
setwd(QUTA_Clust3_wd)
# Read in genind file
QUTA_Clust3_gen <- read.genepop("populations.snps.gen")
# Convert genind file to data.frame
QUTA_Clust3_df <- genind2df(QUTA_Clust3_gen)
# Convert data.frame to genlight object
QUTA_Clust3_genLight <- as.genlight(QUTA_Clust3_df)
# Build a genetic distance matrix from the genlight object
QUTA_Clust3_distMat <- dist(QUTA_Clust3_genLight)
# Build a NJ tree from the genetic distance matrix
QUTA_Clust3_njTree <- nj(QUTA_Clust3_distMat)
# Rename the phylognetic tree, to make more discernible
QUTA_Clust3_njTree$tip.label <- sub('SYST-MOR-000','',QUTA_Clust3_njTree$tip.label)
QUTA_Clust3_njTree$tip.label <- sub('OAK-MOR-00','',QUTA_Clust3_njTree$tip.label)
QUTA_Clust3_njTree$tip.label <- sub('OAK-MOR-000','',QUTA_Clust3_njTree$tip.label)
QUTA_Clust3_njTree$tip.label <- sub('OAK-UMN-000','',QUTA_Clust3_njTree$tip.label)
# Plot
par(mar=c(1,1,2,1)+0.1)
plot(QUTA_Clust3_njTree, 'phylo', main='Clust3 NJ Tree: Abbreviated Sample Names')
# Read in a CSV listing sample names in the order indicated in QUTA_Clust3_njTree$tip.label
clustNames <- read.csv2('/home/akoontz/Documents/QUTA_TxRedOaks/Code/clustIDs_withQUTA_numOrdered.csv', 
                        header = TRUE, sep = ",")[,2]
# Rename and plot
QUTA_Clust3_njTree$tip.label <- clustNames
par(mar=c(1,1,2,1)+0.1)
plot(QUTA_Clust3_njTree, 'phylo', main='Clust3 NJ Tree: K4 Cluster Identities')

# SCY-HYPO-CAN-TAR %%%% ----
QUTA_ScyHypoCanTar_wd <- '/RAID1/QUTA_TX_RedOaks/Genotyping/SNP_Calling/Output/ScyHypoCanTar_seOnly/pop_R97_WL10K/'
setwd(QUTA_ScyHypoCanTar_wd)
# Read in genind file
QUTA_ScyHypoCanTar_gen <- read.genepop("populations.snps.gen")
# Convert genind file to data.frame
QUTA_ScyHypoCanTar_df <- genind2df(QUTA_ScyHypoCanTar_gen)
# Convert data.frame to genlight object
QUTA_ScyHypoCanTar_genLight <- as.genlight(QUTA_ScyHypoCanTar_df)
# Build a genetic distance matrix from the genlight object
QUTA_ScyHypoCanTar_distMat <- dist(QUTA_ScyHypoCanTar_genLight)
# Build a NJ tree from the genetic distance matrix
QUTA_ScyHypoCanTar_njTree <- nj(QUTA_ScyHypoCanTar_distMat)
# Rename the phylognetic tree, to make more discernible
QUTA_ScyHypoCanTar_njTree$tip.label <- sub('SYST-MOR-000','',QUTA_ScyHypoCanTar_njTree$tip.label)
QUTA_ScyHypoCanTar_njTree$tip.label <- sub('OAK-MOR-000','',QUTA_ScyHypoCanTar_njTree$tip.label)
QUTA_ScyHypoCanTar_njTree$tip.label <- sub('OAK-MOR-00','',QUTA_ScyHypoCanTar_njTree$tip.label)
QUTA_ScyHypoCanTar_njTree$tip.label <- sub('OAK-UMN-000','',QUTA_ScyHypoCanTar_njTree$tip.label)
# Plot
par(mar=c(1,1,2,1)+0.1)
plot(QUTA_ScyHypoCanTar_njTree, 'phylo', main='ScyHypoCanTar NJ Tree: Abbreviated Sample Names')

# CLUSTER IDENTITIES
# Make a copy of the original tree
QUTA_ScyHypoCanTar_K6clust_njTree <- QUTA_ScyHypoCanTar_njTree
# Read in a TSV listing sample names/cluster IDs in the order indicated in QUTA_ScyHypoCanTar_njTree$tip.label
clustNames <- 
  read.table('/home/akoontz/Documents/QUTA_TxRedOaks/Code/B_Clustering/3_Clust3/ScyHypoCanTar_SampleList_K6clustIDs.tsv',
             sep = '\t', header = TRUE)
# Replace sample names with abbreviated sample names
clustNames[,1] <- QUTA_ScyHypoCanTar_njTree$tip.label
# Rename and plot
QUTA_ScyHypoCanTar_K6clust_njTree$tip.label <- clustNames[,2]
QUTA_ScyHypoCanTar_K6clust_njTree$tip.label <- paste0(clustNames[,1],'; ',clustNames[,2])
par(mar=c(1,1,2,1)+0.1)
plot(QUTA_ScyHypoCanTar_K6clust_njTree, 'phylo', main='ScyHypoCanTar NJ Tree: K6 Cluster Identities')

# DINSTINCT PHYLO ENTITIES SUBSET %%%% ----
QUTA_SNaQ1_wd <- '/RAID1/QUTA_TX_RedOaks/Genotyping/SNP_Calling/Output/Clust3_seOnly/SNaQ1/'
setwd(QUTA_SNaQ1_wd)
# Read in genind file
QUTA_SNaQ1_gen <- read.genepop("populations.snps.gen")
# Convert genind file to dataframe
QUTA_SNaQ1_df <- genind2df(QUTA_SNaQ1_gen)
# Convert dataframe to genlight object
QUTA_SNaQ1_genLight <- as.genlight(QUTA_SNaQ1_df)
# Build a genetic distance matrix from the genlight object
QUTA_SNaQ1_distMat <- dist(QUTA_SNaQ1_genLight)
# Build a NJ tree from the genetic distance matrix
QUTA_SNaQ1_njTree <- nj(QUTA_SNaQ1_distMat)
# Rename the phylognetic tree, to make more discernible
QUTA_SNaQ1_njTree$tip.label <- sub('SYST-MOR-000','',QUTA_SNaQ1_njTree$tip.label)
QUTA_SNaQ1_njTree$tip.label <- sub('OAK-MOR-000','',QUTA_SNaQ1_njTree$tip.label)
QUTA_SNaQ1_njTree$tip.label <- sub('OAK-MOR-00','',QUTA_SNaQ1_njTree$tip.label)
QUTA_SNaQ1_njTree$tip.label <- sub('OAK-UMN-000','',QUTA_SNaQ1_njTree$tip.label)
# Plot
par(mar=c(1,1,2,1)+0.1)
plot(QUTA_SNaQ1_njTree, 'phylo', main='SNaQ1 NJ Tree: Abbreviated Sample Names')

# SNaQ: ALL TIPS SAMPLE SET %%%% ----
QUTA_SNaQ_AllTips_wd <- '/RAID1/QUTA_TX_RedOaks/Genotyping/SNP_Calling/Output/Clust3_seOnly/SNaQ_AllTips/'
setwd(QUTA_SNaQ_AllTips_wd)
# Read in genind file
QUTA_SNaQ_AllTips_gen <- read.genepop("populations.snps.gen")
# Convert genind file to dataframe
QUTA_SNaQ_AllTips_df <- genind2df(QUTA_SNaQ_AllTips_gen)
# Convert dataframe to genlight object
QUTA_SNaQ_AllTips_genLight <- as.genlight(QUTA_SNaQ_AllTips_df)
# Build a genetic distance matrix from the genlight object
QUTA_SNaQ_AllTips_distMat <- dist(QUTA_SNaQ_AllTips_genLight)
# Build a NJ tree from the genetic distance matrix
QUTA_SNaQ_AllTips_njTree <- nj(QUTA_SNaQ_AllTips_distMat)
# Rename the phylognetic tree, to make more discernible for plot
QUTA_SNaQ_AllTips_njTree$tip.label <- c('EMOR_0361','HYPO_0410','HYPO_1178','MIQU_6235','EMOR_6549',
                                        'GRAV_6557','TARD_6698','TARD_6699','GRAV_7267','TARD_7269',
                                        'TARD_7270','TARD_7271','TARD_7272','TARD_7278','MIQU_7329')
# Plot
par(mar=c(1,1,2,1)+0.1)
plot(QUTA_SNaQ_AllTips_njTree, 'phylo', main='SNaQ_AllTips NJ Tree: Abbreviated Sample Names')
# SNaQ requires an initial topology, to build a phylogenetic network. We write this NJ tree as a Newick
# file, to provide that topology
write.tree(QUTA_SNaQ_AllTips_njTree, 
           file='/RAID1/QUTA_TX_RedOaks/HybridAnalyses/SNaQ/AllTips1/Input/Clust3_SNPs2CF_AllTips.tre')

