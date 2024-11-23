# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% CALCULATE HETEROZYGOSITY %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# This script is used to calculate heterozygosity (Hz) in 3 different groups of species 
# for the Q. tardifolia Texas red oak project. Specifically, it reports this metric 
# in 2 versions of the GravHypoTar dataset, which contains individuals belonging 
# to the "gravesii", "hypoleucoides", and "tardifolia" clusters (although note that 
# individuals with dets other than these cluster names are included within this sample set).
# At the end, t-tests are run to determine whether Hz values are statistically higher in the
# TARD individuals.

library(dartR)
library(report)

# Specify the filepath to the Stacks populations directory
GravHypoTar_Dir <- 
  '/RAID1/QUTA_TX_RedOaks/Genotyping/SNP_Calling/Output/Clust3_seOnly/GravHypoTar_R98'
setwd(GravHypoTar_Dir)

# COMPLETE DATASET, T-TEST ----
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
Hz_Pop <- gl.report.heterozygosity(GravHypoTar_Complete, method = 'pop')

# Report heterozygosity and Fis values by individual
Hz_Ind <- gl.report.heterozygosity(GravHypoTar_Complete, method = 'ind')
# Add a column for cluster IDs, and subset to just Ho values (unclear what f.hom.ref and
# f.hom.alt values are...variances?)
Hz_Ind <- cbind.data.frame(Hz_Ind[,1], as.character(pop(GravHypoTar_Complete)), Hz_Ind[,2])
colnames(Hz_Ind) <- c('SampleNames', 'Clust', 'Ho')
# Write heterozygosities to CSV
write.table(Hz_Ind, 'HzValues_inds.csv', sep=',')
# Extract vectors of heterozygosity values for each cluster ID
TARD_Hz <- Hz_Ind[which(Hz_Ind$Clust=='TARD'),3]
GRAV_Hz <- Hz_Ind[which(Hz_Ind$Clust=='GRAV'),3]
HYPO_Hz <- Hz_Ind[which(Hz_Ind$Clust=='HYPO'),3]

# Perform a t-test for whether TARD Hz values are statistically greater than GRAV's
tTest_TARDvGRAV <- t.test(x=TARD_Hz, y=GRAV_Hz, alternative = 'greater')
# Because p-value is low, we reject the null hypothesis that the difference in means in the
# two groups is 0 (alternative hypothesis is mean is greater in TARD)
report(tTest_TARDvGRAV)
# Perform a t-test for whether TARD Hz values are statistically greater than HYPO's
tTest_TARDvHYPO <- t.test(x=TARD_Hz, y=HYPO_Hz, alternative = 'greater')
# Because p-value is low, we reject the null hypothesis that the difference in means in the
# two groups is 0 (alternative hypothesis is mean is greater in TARD)
report(tTest_TARDvHYPO)

# TRIANGLE PLOTS ----
# Building triangle plots, which are plots of interclass heterozygosity versus hybrid indices,
# and are meant to help distinguish admixture patterns due to hybridization from isolation by distance.
library(triangulaR)
library(vcfR)

# Read in the VCF file generated for 
GHT_vcf <- read.vcfR('populations.snps.vcf')
# Read in a popmap, which is a data.frame with two columns: 'id' and 'pop'
GHT_popmap <- read.table('GravHypoTar_popList.csv', header = TRUE, sep=',')

# Create a new vcfR object composed only of sites above the given allele frequency difference threshold
GHT_vcf_diff <- 
  alleleFreqDiff(vcfR=GHT_vcf, pm=GHT_popmap, p1="GRAV", p2="HYPO", difference=0.9)
# Calculate hybrid index and heterozygosity for each sample. Values are returned in a data.frame
GHT_hybridIndex <- hybridIndex(vcfR=GHT_vcf_diff, pm=GHT_popmap, p1="GRAV", p2="HYPO")

# PLOTTING
# Specify a vector of colors
GHT_cols <- c("#af8dc3", "#7fbf7b", "#bababa", "#878787", "#762a83", "#1b7837")
# Generate a triangle plot
triangle.plot(GHT_hybridIndex, colors=GHT_cols)
