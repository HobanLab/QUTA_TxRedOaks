#!/usr/bin/julia

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% BUILDING PHYLOGENETIC NETWORK OF DISTNCT ENTITIES %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# This is a Julia script used for building a phylogenetic network using SNaQ. This network is used to test
# the likelihood of different phylogenetic placements of Q. tardifolia. The inputs for this script are a 
# table of concordance factors (CF table), which was built using the SNPs2CF tool, and a starting topology,
# which was a simple Newick tree outlining the relations between the samples being analyzed.

# Read in necessary packages
using PhyloNetworks
using PhyloPlots
using Distributed
# Generate multiple nodes, and load PhyloNetworks onto those nodes
addprocs(28)
@everywhere using PhyloNetworks 

# Specify filepath to table of concordance factors, and read in
AllTips_Net1_CF_filepath = "/RAID1/QUTA_TX_RedOaks/HybridAnalyses/SNaQ/AllTips/Input/Clust3_SNPs2CF_AllTips_CFtable.csv"
AllTips_Net1_CF = readTableCF(AllTips_Net1_CF_filepath)
# Specify filepath to the starting topology (a Newick text file), and read in
AllTips_Net1_Tree_filepath = "/RAID1/QUTA_TX_RedOaks/HybridAnalyses/SNaQ/AllTips/Input/Clust3_SNPs2CF_AllTips.tre"
AllTips_Net1_Tree = readTopology(AllTips_Net1_Tree_filepath)
# Specify output directory
AllTips_Out = "/RAID1/QUTA_TX_RedOaks/HybridAnalyses/SNaQ/AllTips/Output/AllTips_Net1"
# Call SNaQ to build the initial phylogenetic network. The hmax argument is set to 1 
# The default for SNaQ is to use 10 runs.
AllTips_Net1 = snaq!(AllTips_Net1_Tree,  AllTips_Net1_CF, hmax=1, filename=AllTips_Out, seed=123)

# Plot the resulting phylogenetic network and save to a PDF file
using RCall
R"pdf"("AllTips_Net1.pdf", width=15, height=9);
plot(AllTips_Net1, :R);
R"dev.off()";

?snaq!