# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% CALCULATE HETEROZYGOSITY %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

library(adegenet)

# Read popgen file and calculate heterozygosity
GravHypoTar.gen <- read.genepop('/RAID1/QUTA_TX_RedOaks/Genotyping/SNP_Calling/Output/Clust3_seOnly/GravHypoTar_R98_WL10K/populations.snps.gen')

pop(GravHypoTar.gen)

indNames(GravHypoTar.gen)

Hs()
