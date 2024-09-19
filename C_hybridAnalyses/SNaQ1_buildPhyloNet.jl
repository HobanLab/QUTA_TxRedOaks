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

# Read in the table of concordance factors
Clust3_CF_8sp = readTableCF("Input/Clust3_SNaQ1_CFtable_8sp.csv")
# Specify the starting topology (a Newick text file)
Clust3_Tree = readTopology("Input/SNaQ1_PhyloEnts.tre")
# Call SNaQ to build the initial phylogenetic network. The hmax argument is set to 0 because we're initially 
# exploring no hybridization events. The default for SNaQ is to use 10 runs.
Net0 = snaq!(Clust3_Tree,  Clust3_CF_8sp, hmax=0, filename="Output/Net0_SNaQ1_8sp", seed=123)

