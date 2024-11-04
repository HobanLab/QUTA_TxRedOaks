# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% CALCULATE HETEROZYGOSITY %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# This script is used to calculate heterozygosity in 3 different
# groups of species for the Q. tardifolia Texas red oak project.
# Specifically, it reports this metric in 2 versions of the 
# GravHypoTar dataset, which contains individuals belonging to the 
# "gravesii", "hypoleucoides", and "tardifolia" clusters (although
# note that individuals with dets other than these cluster names
# are included within this sample set).

# 1. The "complete" dataset, using all of the loci present at R98
# 2. The "WL10K" dataset, using 10,000 randomly whitelisted loci

library(dartR)

# COMPLETE DATASET ----
# Specify the filepath to the Stacks populations directory
GravHypoTar_Dir <- 
  '/RAID1/QUTA_TX_RedOaks/Genotyping/SNP_Calling/Output/Clust3_seOnly/GravHypoTar_R98'
setwd(GravHypoTar_Dir)

# Specify the filepath to the relevant VCF file, and read it in
GravHypoTar_Complete <- gl.read.vcf('populations.snps.vcf')
# Specify the populations of the samples. This is necessary because
# population assignments aren't stored within the VCF. Pop IDs are
# first read in from a CSV to a slot in the genlight object
GravHypoTar_Complete$other$ind.metrics <- 
  read.table('GravHypoTar_popList.csv', header=T, sep = ',')
pop(GravHypoTar_Complete) <- GravHypoTar_Complete$other$ind.metrics$pop
# Check that population assignments look good
cbind(indNames(GravHypoTar_Complete), as.character(pop(GravHypoTar_Complete)))
# Report heterozygosity and Fis values by population
gl.report.heterozygosity(GravHypoTar_Complete, method = 'pop')
# Report heterozygosity and Fis values by individual
gl.report.heterozygosity(GravHypoTar_Complete, method = 'ind')
