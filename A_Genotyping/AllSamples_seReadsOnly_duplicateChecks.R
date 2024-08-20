# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% QUTA & TEXAS RED OAKS : ANALYZING GENEPOP FILES %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# This is an exploratory script, for checking out the number of loci that are distinct between replicate samples

library(adegenet)
QUTA_TRO_wd <- '/RAID1/QUTA_TX_RedOaks/Genotyping/SNP_Calling/Output/All_seReadsOnly/'

# %%% R97 DATASET ----
# Read in genepop file for All samples (147), single-end data only, R97
R97_genind <- read.genepop(paste0(QUTA_TRO_wd,"pop_R97/populations.snps.gen"))
# Duplicates: SYST-MOR-0007278 and SYST-MOR-0007278rpt
R97_genind@tab[156:157, 1:5]
# Which alleles are different between replicates
which(R97_genind@tab[156,] != R97_genind@tab[157,])
# How many alleles are different between replicates (54)
length(which(R97_genind@tab[156,] != R97_genind@tab[157,]))
# What proportion of alleles are different between replicates (54/24,246 = 0.22%)
length(which(R97_genind@tab[156,] != R97_genind@tab[157,]))/ncol(R97_genind@tab)

# %%% R98 DATASET ----
# Read in genepop file for All samples (147), single-end data only, R98
R98_genind <- read.genepop(paste0(QUTA_TRO_wd,"pop_R98/populations.snps.gen"))
# Duplicates: SYST-MOR-0007278 and SYST-MOR-0007278rpt
R98_genind@tab[156:157, 1:5]
# Which alleles are different between replicates
which(R98_genind@tab[156,] != R98_genind@tab[157,])
# How many alleles are different between replicates (12)
length(which(R98_genind@tab[156,] != R98_genind@tab[157,]))
# What proportion of alleles are different between replicates (12/4,030 = 0.30%)
length(which(R98_genind@tab[156,] != R98_genind@tab[157,]))/ncol(R98_genind@tab)
