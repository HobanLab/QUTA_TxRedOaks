# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% QUERCUS TARDIFOLIA F3 STATISTICS %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# This script is used to calculate F3 statistics for the set of individuals utilized in Clustering run #3
# of the QUTA Texas red oak population analyses. It reads in the PLINK file from that run, and then
# generates a table of F3 statistics using the admixtools library

library(admixtools)

# %%% CALCULATING F3 STATS FOR TARDIFOLIA %%% ----
# All populations 'wild'
plinkDir <- '/RAID1/QUTA_TX_RedOaks/Genotyping/SNP_Calling/Output/Clust3_seOnly/pop_R98/plinkFiles/'
setwd(plinkDir)

Clust3_seOnly_R98_f2Arr <- f2_from_geno(pref='Clust3_seOnly_R98', verbose = TRUE)
# RESULT: Error in extract_f2(pref, outdir = NULL, inds = inds, pops = pops, ...) : There are no informative SNPs!
# NOTE: unclear if this is because of the way the populations argument is specified, or if there are simply no
#       SNPs that have 0 missing data (i.e. SNPs present in every single individuals)

# Each sample its own population
plinkDir <- '/RAID1/QUTA_TX_RedOaks/Genotyping/SNP_Calling/Output/Clust3_seOnly/pop_R98_unqPops/plinkFiles/'
setwd(plinkDir)

Clust3_seOnly_R98_f2Arr <- f2_from_geno(pref='Clust3_seOnly_R98', verbose = TRUE)
# NOTE: Warning--no genetic linkage map found.
# NOTE: the third dimension is referred to as "each SNP block". There are 43 in our data, which suggest that 
#       our data is split into 43 SNP blocks

# If pop argument for qp3pop isn't specified, all population combinations will be calculated, which would be...overkill
# Start by providing population combinations for a subset of samples: one tardifolia, one gravesii, one hypoleucoides
QUTA_testPops_A <- c('7271','1178','6557')
Clust3_seOnly_R98_f3dfA <- f3(data=Clust3_seOnly_R98_f2Arr, pop1 = QUTA_testPops_A)

# Now, perhaps we can explore using a list of all tardifolia, gravesii, and hypoleucoides individuals?
QUTA_tardifolia_Pops <- c('6698','6699','7269','7270','7271','7272','7278')
QUTA_gravesii_Pops <- c('6628','6553','6556','6557','6558','7262','7267')
QUTA_hypoleucoides_Pops <- c('6636','7274','7275','410','416','180','896')
QUTA_emoryi_Pops <- c('6533','6550')

# Calculate F3 stats for QUTA and expected parentals (gravesii and hypoleucoides)
Clust3_seOnly_R98_f3df_QUTA_Exp <- f3(data=Clust3_seOnly_R98_f2Arr, pop1 = QUTA_tardifolia_Pops, 
                                      pop2=QUTA_gravesii_Pops, pop3=QUTA_hypoleucoides_Pops)
# Write tibble results to a table
write.table(Clust3_seOnly_R98_f3df_QUTA_Exp, 
            file = '/home/akoontz/Documents/QUTA_TxRedOaks/Code/QUTA_F3Results_ExpectedParentals.csv', sep = ',')

# Calculate F3 stats for QUTA and all possible parentals (as a test)
Clust3_seOnly_R98_f3df_QUTA <- f3(data=Clust3_seOnly_R98_f2Arr, pop1=QUTA_tardifolia_Pops, 
                                  pop2=colnames(Clust3_seOnly_R98_f2Arr), pop3=colnames(Clust3_seOnly_R98_f2Arr))
# Write tibble results to a table
write.table(Clust3_seOnly_R98_f3df_QUTA, 
            file = '/home/akoontz/Documents/QUTA_TxRedOaks/Code/QUTA_F3Results_AllCombinations.csv', sep = ',')

# %%% TESTING HYBRIDIZATION IN EMORYI %%% ----
Clust3_seOnly_R98_f3df_emoryi <- f3(data=Clust3_seOnly_R98_f2Arr, pop1 = QUTA_emoryi_Pops, 
                                    pop2=QUTA_gravesii_Pops, pop3=QUTA_hypoleucoides_Pops)
print(n=98, Clust3_seOnly_R98_f3df_emoryi)
# %%% TESTING HYBRIDIZATION IN CANBYI %%% ----
f3(data=Clust3_seOnly_R98_f2Arr, pop1='7351', pop2='7328', pop3 ='7267')
# DO THIS MORE: Look into support for hybridization patterns in canbyi, others

# %%% CALCULATE F3 VALUES FOR ALL POSSIBLE COMBINATIONS %%% ----
# This command actually doesn't take too long to run--less than 30 seconds
Clust3_seOnly_R98_f3df_Total <- f3(data=Clust3_seOnly_R98_f2Arr)

# %%% EXPLORING F4 VALUES %%% ----
# Low f4 values (when samples A and B come from a totally different, not admixed cluster than C and D)? Sort of...
f4(data=Clust3_seOnly_R98_f2Arr, pop1 = '6628', pop2='6553', pop3='7329', pop4='7328')
f4(data=Clust3_seOnly_R98_f2Arr, pop1 = '6628', pop2='6553', pop3='6636', pop4='7274')
# High f4 values (when all samples come from the same cluster)? No
f4(data=Clust3_seOnly_R98_f2Arr, pop1 = '6628', pop2='6553', pop3='6556', pop4='6557')
f4(data=Clust3_seOnly_R98_f2Arr, pop1 = '410', pop2='416', pop3='180', pop4='896')
# So, it's difficult to interpret the magnitude of 4 values...

# %%% EXPLORING ADMIXTURE GRAPHS %%% ----
find_graphs(Clust3_seOnly_R98_f2Arr, numadmix = 3, verbose = TRUE)
