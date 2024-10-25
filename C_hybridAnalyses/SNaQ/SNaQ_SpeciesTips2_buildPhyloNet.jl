#!/usr/bin/julia

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% BUILDING PHYLOGENETIC NETWORK OF DISTNCT ENTITIES %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# This is a Julia script used for building a phylogenetic network using SNaQ. This network is used to test
# the likelihood of different phylogenetic placements of Q. tardifolia. The inputs for this script are a 
# table of concordance factors (CF table), which was built using the SNPs2CF tool, and a starting topology,
# which was a simple Newick tree outlining the relations between the 5 species being analyzed.

# This script also takes the step of mapping the CF table such that individuals are grouped
# into species. This is done prior to calling the readTableCF command.

# Read in necessary packages
using PhyloNetworks
using PhyloPlots
using Distributed
using CSV, DataFrames
# Generate multiple nodes, and load PhyloNetworks onto those nodes
addprocs(16)
@everywhere using PhyloNetworks 
# Specify the filepath containing the SNaQ input files
SNaQ_path = "/RAID1/QUTA_TX_RedOaks/HybridAnalyses/SNaQ/SpeciesTips/Input/"

# Read in the mapping file and the CF table
CFmappingFile_filepath = joinpath(SNaQ_path,"SNaQ_sampleMap.csv")
SpeciesTips_CF_filepath = joinpath(SNaQ_path,"Clust3_SNPs2CF_SpeciesTips_CFtable.csv")
# Merge individuals into species, using the mapping file
New_SpeciesTips_CF = mapAllelesCFtable(CFmappingFile_filepath, SpeciesTips_CF_filepath)
SpeciesTips_CF = readTableCF!(New_SpeciesTips_CF, mergerows=true)

# Specify filepath to the starting topology (a Newick text file), and read in
SpeciesTips_Tree_filepath = "/RAID1/QUTA_TX_RedOaks/HybridAnalyses/SNaQ/SpeciesTips/Input/Clust3_SpeciesTips.tre"
SpeciesTips_Tree = readTopology(SpeciesTips_Tree_filepath)
# Specify output directory
SpeciesTips_Net0_Out = "/RAID1/QUTA_TX_RedOaks/HybridAnalyses/SNaQ/SpeciesTips/Output/SpeciesTips_Net0"
SpeciesTips_Net1_Out = "/RAID1/QUTA_TX_RedOaks/HybridAnalyses/SNaQ/SpeciesTips/Output/SpeciesTips_Net1"
SpeciesTips_Net2_Out = "/RAID1/QUTA_TX_RedOaks/HybridAnalyses/SNaQ/SpeciesTips/Output/SpeciesTips_Net2"

# ---- RUNNING SNAQ ----
# Call SNaQ to build the initial phylogenetic network, varying the hmax value. 
# The default for SNaQ is to use 10 runs.
# NET 0
SpeciesTips_Net0 = snaq!(SpeciesTips_Tree,  SpeciesTips_CF, hmax=0, filename=SpeciesTips_Net0_Out, seed=123)
# Plot the resulting phylogenetic network and save to a PDF file
using RCall
R"pdf"("SpeciesTips_Net0.pdf", width=15, height=9);
plot(SpeciesTips_Net0, :R);
R"dev.off()";

# NET 1
SpeciesTips_Net1 = snaq!(SpeciesTips_Tree,  SpeciesTips_CF, hmax=1, filename=SpeciesTips_Net1_Out, seed=456)
# Plot the resulting phylogenetic network and save to a PDF file
using RCall
R"pdf"("SpeciesTips_Net1.pdf", width=15, height=9);
plot(SpeciesTips_Net1, :R);
R"dev.off()";

# NET 2
SpeciesTips_Net2 = snaq!(SpeciesTips_Tree,  SpeciesTips_CF, hmax=2, filename=SpeciesTips_Net2_Out, seed=789)
# Plot the resulting phylogenetic network and save to a PDF file
using RCall
R"pdf"("SpeciesTips_Net2.pdf", width=15, height=9);
plot(SpeciesTips_Net2, :R);
R"dev.off()";
