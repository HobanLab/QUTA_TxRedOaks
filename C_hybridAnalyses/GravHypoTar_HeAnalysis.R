# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% CALCULATE HETEROZYGOSITY %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# This script is used to calculate expected heterozygosity in 3 different
# groups of species for the Q. tardifolia Texas red oak project.
# Specifically, it reports this metric in 2 versions of the 
# GravHypoTar dataset, which contains individuals belonging to the 
# "gravesii", "hypoleucoides", and "tardifolia" clusters (although
# note that individuals with dets other than these cluster names
# are included within this sample set).

# 1. The "complete" dataset, using all of the loci present at R98
# 2. The "WL10K" dataset, using 10,000 randomly whitelisted loci

library(adegenet)
# Specify the filepath to the Stacks populations directory
GravHypoTar_Dir <- 
  '/RAID1/QUTA_TX_RedOaks/Genotyping/SNP_Calling/Output/Clust3_seOnly/'

# COMPLETE DATASET ----
# Specify the filepath to the relevant genepop file, and read it in
GravHypoTar_Complete_filepath <- 
  paste0(GravHypoTar_Dir, 'GravHypoTar_R98/populations.snps.gen')
GravHypoTar.Complete.gen <- read.genepop(GravHypoTar_Complete_filepath)
# Specify the populations of each sample. These are in such an obscure 
# order because samples in the genepop file are ordered alphabetically
pop(GravHypoTar.Complete.gen) <- 
  c(rep('HYPO',8), rep('GRAV',12),'HYPO', rep('TARD',2), 'HYPO',
    rep('TARD',2), 'HYPO', rep('TARD',4), 'GRAV', 'HYPO', rep('GRAV',2),
    rep('TARD',4), rep('HYPO',2), 'TARD', 'HYPO')
# Calculate heterosygosity
Hs(GravHypoTar.Complete.gen)

# WL10K DATASET ----
# Specify the filepath to the relevant genepop file, and read it in
GravHypoTar_WL10K_filepath <- 
  paste0(GravHypoTar_Dir, 'GravHypoTar_R98_WL10K/populations.snps.gen')
GravHypoTar.WLK10K.gen <- read.genepop(GravHypoTar_WL10K_filepath)
# Specify the populations of each sample. These are in such an obscure 
# order because samples in the genepop file are ordered alphabetically
pop(GravHypoTar.WLK10K.gen) <- 
  c(rep('HYPO',8), rep('GRAV',12),'HYPO', rep('TARD',2), 'HYPO',
    rep('TARD',2), 'HYPO', rep('TARD',4), 'GRAV', 'HYPO', rep('GRAV',2),
    rep('TARD',4), rep('HYPO',2), 'TARD', 'HYPO')
# Calculate heterosygosity
Hs(GravHypoTar.WLK10K.gen)
