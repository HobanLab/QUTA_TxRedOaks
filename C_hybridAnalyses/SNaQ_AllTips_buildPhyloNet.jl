#!/usr/bin/julia

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% BUILDING PHYLOGENETIC NETWORK OF DISTNCT ENTITIES %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# This is a Julia script used for building a phylogenetic network using SNaQ. This network is used to test
# the likelihood of different phylogenetic placements of Q. tardifolia. The inputs for this script are a 
# table of concordance factors (CF table), which was built using the SNPs2CF tool, and a starting topology,
# which was a simple Newick tree outlining the relations between the 8 samples being analyzed.

# Read in necessary packages
using PhyloNetworks
using PhyloPlots
using Distributed
# Generate multiple nodes, and load PhyloNetworks onto those nodes
addprocs(28)
@everywhere using PhyloNetworks 

# Specify filepath to table of concordance factors, and read in
SNaQ_AllTips_15sp_CF_filepath = "/RAID1/QUTA_TX_RedOaks/HybridAnalyses/SNaQ/AllTips1/Input/Clust3_SNPs2CF_AllTips_CFtable.csv"
SNaQ_AllTips_15sp_CF = readTableCF(SNaQ_AllTips_15sp_CF_filepath)
# Specify filepath to the starting topology (a Newick text file), and read in
SNaQ_AllTips_15sp_Tree_filepath = "/RAID1/QUTA_TX_RedOaks/HybridAnalyses/SNaQ/AllTips1/Input/Clust3_SNPs2CF_AllTips.tre"
SNaQ_AllTips_15sp_Tree = readTopology(SNaQ_AllTips_15sp_Tree_filepath)
# Specify output directory
SNaQ_AllTips_Out = "/RAID1/QUTA_TX_RedOaks/HybridAnalyses/SNaQ/AllTips1/Output/Net1_AllTips_15sp"
# Call SNaQ to build the initial phylogenetic network. The hmax argument is set to 0 because we're initially 
# exploring no hybridization events. The default for SNaQ is to use 10 runs.
SNaQ_AllTips_Net1 = snaq!(SNaQ_AllTips_15sp_Tree,  SNaQ_AllTips_15sp_CF, hmax=1, filename=SNaQ_AllTips_Out, seed=123)

# Plot the resulting phylogenetic network and save to a PDF file
using RCall
R"pdf"("Plot-AllTips_Net1.pdf", width=5, height=3);
plot(SNaQ_AllTips_Net1, :R);
R"dev.off()";
