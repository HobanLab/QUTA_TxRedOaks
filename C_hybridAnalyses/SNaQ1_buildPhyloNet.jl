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
addprocs(10)
@everywhere using PhyloNetworks 

# Specify filepath to table of concordance factors, and read in
Clust3_CF_8sp_filepath = "/RAID1/QUTA_TX_RedOaks/HybridAnalyses/SNaQ/SNaQ1/Input/Clust3_SNaQ1_CFtable_8sp.csv"
Clust3_CF_8sp = readTableCF(Clust3_CF_8sp_filepath)
# Specify filepath to the starting topology (a Newick text file), and read in
Clust3_Tree_filepath = "/RAID1/QUTA_TX_RedOaks/HybridAnalyses/SNaQ/SNaQ1/Input/SNaQ1_PhyloEnts.tre"
Clust3_Tree = readTopology(Clust3_Tree_filepath)
# Specify output directory
SNaQ1_Out = "/RAID1/QUTA_TX_RedOaks/HybridAnalyses/SNaQ/SNaQ1/Output/Net0_SNaQ1_8sp"
# Call SNaQ to build the initial phylogenetic network. The hmax argument is set to 0 because we're initially 
# exploring no hybridization events. The default for SNaQ is to use 10 runs.
SNaQ1_Net0 = snaq!(Clust3_Tree,  Clust3_CF_8sp, hmax=0, filename=SNaQ1_Out, seed=123)

# Plot the resulting phylogenetic network and save to a PDF file
using RCall
R"pdf"("Plot-Net0.pdf", width=5, height=3);
plot(SNaQ1_Net0, :R);
R"dev.off()";
