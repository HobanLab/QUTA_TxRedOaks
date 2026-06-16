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

library(dartR, quietly = TRUE)
library(report)

# %%%% ROUND 2, OCTOBER 2025: DOWNSAMPLED DATASET %%%% ---- 
# Specify the filepath to the Stacks populations directory for the Clust 3 dataset
# The outputs in this directory utilize a popmap file that includes 3 groups:
# 7 gravesii, 7 hypoleucoides, and 7 tardifolia
GravHypoTar_Dir <- 
  '/RAID1/QUTA_TX_RedOaks/Genotyping/SNP_Calling/Output/Clust3_seOnly/GravHypoTar2_R98'
setwd(GravHypoTar_Dir)
# Variable file path: directory to save images to
outputPlotDir <- 
  "/home/akoontz/Documents/QUTA_TxRedOaks/Documentation/Images/ManuscriptDraft_IJPS/"

# T-TEST ----
# Specify the filepath to the relevant VCF file, and read it in
GravHypoTar_Complete <- gl.read.vcf('populations.snps.vcf')
# Specify the populations of the samples. This is necessary because
# population assignments aren't stored within the VCF. Pop IDs are
# first read in from a CSV to a slot in the genlight object
GravHypoTar_Complete$other$ind.metrics <- 
  read.table('GravHypoTar2_popList.csv', header=T, sep = ',')
pop(GravHypoTar_Complete) <- GravHypoTar_Complete$other$ind.metrics$pop
# Check that population assignments look good
cbind(indNames(GravHypoTar_Complete), as.character(GravHypoTar_Complete@pop))
# Report heterozygosity and Fis values by population
Hz_Pop <- gl.report.heterozygosity(GravHypoTar_Complete, method = 'pop')

# Report heterozygosity and Fis values by individual
Hz_Ind <- gl.report.heterozygosity(GravHypoTar_Complete, method = 'ind')
# Add a column for cluster IDs, and subset to just Ho values (unclear what f.hom.ref and
# f.hom.alt values are...variances?)
Hz_Ind <- cbind.data.frame(Hz_Ind[,1], as.character(GravHypoTar_Complete@pop), Hz_Ind[,2])
colnames(Hz_Ind) <- c('SampleNames', 'Clust', 'Ho')
# Write heterozygosities to CSV
write.table(Hz_Ind, 'HzValues_Downsampled_inds.csv', sep=',')
# Extract vectors of heterozygosity values for each cluster ID
TARD_Hz <- Hz_Ind[which(Hz_Ind$Clust=='TARD'),3]
GRAV_Hz <- Hz_Ind[which(Hz_Ind$Clust=='GRAV'),3]
HYPO_Hz <- Hz_Ind[which(Hz_Ind$Clust=='HYPO'),3]

# Perform a t-test for whether TARD Hz values are statistically greater than GRAV's
# Argument var.equal=FALSE (default) specifies Welch's t-test (doesn't assume equal variances)
tTest_TARDvGRAV <- t.test(x=TARD_Hz, y=GRAV_Hz, alternative = 'greater')
# Because p-value is low, we reject the null hypothesis that the difference in means in the
# two groups is 0 (alternative hypothesis is mean is greater in TARD)
report(tTest_TARDvGRAV)
# Perform a t-test for whether TARD Hz values are statistically greater than HYPO's
# Argument var.equal=FALSE (default) specifies Welch's t-test (doesn't assume equal variances)
tTest_TARDvHYPO <- t.test(x=TARD_Hz, y=HYPO_Hz, alternative = 'greater')
# Because p-value is low, we reject the null hypothesis that the difference in means in the
# two groups is 0 (alternative hypothesis is mean is greater in TARD)
report(tTest_TARDvHYPO)

# TRIANGLE PLOTS ----
# Building triangle plots, which are plots of interclass heterozygosity versus hybrid indices,
# and are meant to help distinguish admixture patterns due to hybridization from isolation by distance.
library(triangulaR)
library(vcfR)

# Read in the VCF file generated for heterozygosity comparisons
GHT_vcf <- read.vcfR('populations.snps.vcf')
# Read in a popmap, which is a data.frame with two columns: 'id' and 'pop'
GHT_popmap <- read.table('GravHypoTar2_popList.csv', header = TRUE, sep=',')
# Specify a value for the allele frequency difference threshold (only alleles with frequency differences
# of this amount or greater in parental groups are considered for calculating interclass heterozygosity)
alleleFreqThresh <- 0.75
# Create a new vcfR object composed only of sites above the given allele frequency difference threshold
GHT_vcf_diff_0.75 <- 
  alleleFreqDiff(vcfR=GHT_vcf, pm=GHT_popmap, p1="GRAV", p2="HYPO", difference=alleleFreqThresh)
# Calculate hybrid index and heterozygosity for each sample. Values are returned in a data.frame
GHT_hybridIndex <- hybridIndex(vcfR=GHT_vcf_diff_0.75, pm=GHT_popmap, p1="GRAV", p2="HYPO")

# PLOTTING
# Specify a vector of colors
GHT_cols <- c("#af8dc3", "#7fbf7b", "#bababa", "#878787", "#762a83", "#1b7837")
# Generate a triangle plot
tri_plot_75 <- triangulaR::triangle.plot(GHT_hybridIndex, colors=GHT_cols, cex=3)
# Update legend on the plot
tri_plot_75 <- tri_plot_75 + theme(legend.position = c(0.85, 0.7)) +
  scale_color_manual(
    name   = "Taxa",
    values = c("GRAV" = GHT_cols[[1]],
               "HYPO" = GHT_cols[[2]],
               "TARD" = GHT_cols[[3]]),
    labels = c("Q. gravesii", "Q. hypoleucoides", "Q. tardifolia")
  )
# Print the allele frequency threshold value and the number of differentiated sites
tri_plot_75 <- tri_plot_75 + 
  annotate('text', x=0.1, y=0.9, label=paste0('Allele frequency threshold: ', alleleFreqThresh)) + 
  annotate('text', x=0.1, y=0.8, label=paste0('Number of sites: ', nrow(GHT_vcf_diff_0.75@fix)))
print(tri_plot_75)
# Save the image to a TIFF
tiff(file = paste0(outputPlotDir, "FigS11_triPlot_75.tiff"), width = 1115, height = 389)
print(tri_plot_75)
dev.off()

# Repeat the process, but using a higher allele frequency threshold (0.9)
alleleFreqThresh <- 0.9
GHT_vcf_diff_0.9 <- 
  alleleFreqDiff(vcfR=GHT_vcf, pm=GHT_popmap, p1="GRAV", p2="HYPO", difference=alleleFreqThresh)
GHT_hybridIndex <- hybridIndex(vcfR=GHT_vcf_diff_0.9, pm=GHT_popmap, p1="GRAV", p2="HYPO")
# Generate a triangle plot
tri_plot_90 <- triangulaR::triangle.plot(GHT_hybridIndex, colors=GHT_cols, cex=3)
# Update legend on the plot
tri_plot_90 <- tri_plot_90 + theme(legend.position = c(0.85, 0.7)) +
  scale_color_manual(
    name   = "Taxa",
    values = c("GRAV" = GHT_cols[[1]],
               "HYPO" = GHT_cols[[2]],
               "TARD" = GHT_cols[[3]]),
    labels = c("Q. gravesii", "Q. hypoleucoides", "Q. tardifolia")
  )
# Print the allele frequency threshold value and the number of differentiated sites
tri_plot_90 <- tri_plot_90 + 
  annotate('text', x=0.1, y=0.9, label=paste0('Allele frequency threshold: ', alleleFreqThresh)) + 
  annotate('text', x=0.1, y=0.8, label=paste0('Number of sites: ', nrow(GHT_vcf_diff_0.9@fix)))
print(tri_plot_90)
# Save the image to a TIFF
tiff(file = paste0(outputPlotDir, "Fig5_triPlot_90.tiff"), width = 1115, height = 389)
print(tri_plot_90)
dev.off()

# FILE CONVERSION STEPS, FOR NEWHYBRIDS ANALYSIS
# The commands below are used to convert the VCF of 211 fixed differences to a
# genepop object. This was used in a separate analysis.
write.vcf(GHT_vcf_diff_0.9, file='NewHybridsInputs/QUTA_GravHypoTar2_R98_AFT09.vcf')
# Create genlight from VCF, and then convert genlight to genepop file
QUTA_GravHypoTar2_AFT0.9_GL <- gl.read.vcf('NewHybridsInputs/QUTA_GravHypoTar2_R98_AFT09.vcf')
QUTA_GravHypoTar2_AFT0.9_GP <- 
  gl2genepop(QUTA_GravHypoTar2_AFT0.9_GL, outfile = 'QUTA_GravHypoTar2_R98_AFT09.gen',
             outpath = './NewHybridsInputs/')

# Repeat the process, but using a higher allele frequency threshold (1.0)
alleleFreqThresh <- 1
GHT_vcf_diff_1 <-
  alleleFreqDiff(vcfR=GHT_vcf, pm=GHT_popmap, p1="GRAV", p2="HYPO", difference=alleleFreqThresh)
GHT_hybridIndex <- hybridIndex(vcfR=GHT_vcf_diff_1, pm=GHT_popmap, p1="GRAV", p2="HYPO")
# Generate a triangle plot
tri_plot_100 <- triangulaR::triangle.plot(GHT_hybridIndex, colors=GHT_cols, cex=3)
# Update legend on the plot
tri_plot_100 <- tri_plot_100 + theme(legend.position = c(0.85, 0.7)) +
  scale_color_manual(
    name   = "Taxa",
    values = c("GRAV" = GHT_cols[[1]],
               "HYPO" = GHT_cols[[2]],
               "TARD" = GHT_cols[[3]]),
    labels = c("Q. gravesii", "Q. hypoleucoides", "Q. tardifolia")
  )
# Print the allele frequency threshold value and the number of differentiated sites
tri_plot_100 <- tri_plot_100 + 
  annotate('text', x=0.1, y=0.9, label=paste0('Allele frequency threshold: ', alleleFreqThresh)) + 
  annotate('text', x=0.1, y=0.8, label=paste0('Number of sites: ', nrow(GHT_vcf_diff_1@fix)))
print(tri_plot_100)
# Save the image to a TIFF
tiff(file = paste0(outputPlotDir, "FigS12_triPlot_100.tiff"), width = 1115, height = 389)
print(tri_plot_100)
dev.off()

# FILE CONVERSION STEPS, FOR NEWHYBRIDS ANALYSIS
# The commands below are used to convert the VCF of 211 fixed differences to a
# genepop object. This was used in a separate analysis.
write.vcf(GHT_vcf_diff_1, file='NewHybridsInputs/QUTA_GravHypoTar2_R98_AFT1.vcf')
# Create genlight from VCF, and then convert genlight to genepop file
QUTA_GravHypoTar2_AFT1_GL <- gl.read.vcf('NewHybridsInputs/QUTA_GravHypoTar2_R98_AFT1.vcf')
QUTA_GravHypoTar2_AFT1_GP <- 
  gl2genepop(QUTA_GravHypoTar2_AFT1_GL, outfile = 'QUTA_GravHypoTar2_R98_AFT1.gen',
             outpath = './NewHybridsInputs/')

# %%%% ROUND 1, NOVEMBER 2024: COMPLETE DATASET %%%% ---- 
# Specify the filepath to the Stacks populations directory for the Clust 3 dataset
# The outputs in this directory utilize a popmap file that includes 3 (general) groups:
# 15 gravesii (some "aff. gravesii", some one "graciliformis (ish)")
# 15 hypoleucoides (and, notably, scytophylla individuals...)
# 13 tardifolia (as well as known gravesii x hypoleucoides hybrids...)
GravHypoTar_Dir <- 
  '/RAID1/QUTA_TX_RedOaks/Genotyping/SNP_Calling/Output/Clust3_seOnly/GravHypoTar_R98'
setwd(GravHypoTar_Dir)
# T-TEST ----
# Specify the filepath to the relevant VCF file, and read it in
GravHypoTar_Complete <- gl.read.vcf('populations.snps.vcf')
# Specify the populations of the samples. This is necessary because
# population assignments aren't stored within the VCF. Pop IDs are
# first read in from a CSV to a slot in the genlight object
GravHypoTar_Complete$other$ind.metrics <- 
  read.table('GravHypoTar_popList.csv', header=T, sep = ',')
pop(GravHypoTar_Complete) <- GravHypoTar_Complete$other$ind.metrics$pop
# Check that population assignments look good
cbind(indNames(GravHypoTar_Complete), as.character(GravHypoTar_Complete@pop))
# Report heterozygosity and Fis values by population
Hz_Pop <- gl.report.heterozygosity(GravHypoTar_Complete, method = 'pop')

# Report heterozygosity and Fis values by individual
Hz_Ind <- gl.report.heterozygosity(GravHypoTar_Complete, method = 'ind')
# Add a column for cluster IDs, and subset to just Ho values (unclear what f.hom.ref and
# f.hom.alt values are...variances?)
Hz_Ind <- cbind.data.frame(Hz_Ind[,1], as.character(GravHypoTar_Complete@pop), Hz_Ind[,2])
colnames(Hz_Ind) <- c('SampleNames', 'Clust', 'Ho')
# Write heterozygosities to CSV
write.table(Hz_Ind, 'HzValues_Complete_inds.csv', sep=',')
# Extract vectors of heterozygosity values for each cluster ID
TARD_Hz <- Hz_Ind[which(Hz_Ind$Clust=='TARD'),3]
GRAV_Hz <- Hz_Ind[which(Hz_Ind$Clust=='GRAV'),3]
HYPO_Hz <- Hz_Ind[which(Hz_Ind$Clust=='HYPO'),3]

# Perform a t-test for whether TARD Hz values are statistically greater than GRAV's
# Argument var.equal=FALSE (default) specifies Welch's t-test (doesn't assume equal variances)
tTest_TARDvGRAV <- t.test(x=TARD_Hz, y=GRAV_Hz, alternative = 'greater')
# Because p-value is low, we reject the null hypothesis that the difference in means in the
# two groups is 0 (alternative hypothesis is mean is greater in TARD)
report(tTest_TARDvGRAV)
# Perform a t-test for whether TARD Hz values are statistically greater than HYPO's
# Argument var.equal=FALSE (default) specifies Welch's t-test (doesn't assume equal variances)
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
# Specify a value for the allele frequency difference threshold (only alleles with frequency differences
# of this amount or greater in parental groups are considered for calculating interclass heterozygosity)
alleleFreqThresh <- 0.75
# Create a new vcfR object composed only of sites above the given allele frequency difference threshold
GHT_vcf_diff <- 
  alleleFreqDiff(vcfR=GHT_vcf, pm=GHT_popmap, p1="GRAV", p2="HYPO", difference=alleleFreqThresh)
# Calculate hybrid index and heterozygosity for each sample. Values are returned in a data.frame
GHT_hybridIndex <- hybridIndex(vcfR=GHT_vcf_diff, pm=GHT_popmap, p1="GRAV", p2="HYPO")

# PLOTTING
# Specify a vector of colors
GHT_cols <- c("#af8dc3", "#7fbf7b", "#bababa", "#878787", "#762a83", "#1b7837")
# Generate a triangle plot
tri_plot <- triangulaR::triangle.plot(GHT_hybridIndex, colors=GHT_cols, cex=3)
# Print the allele frequency threshold value and the number of differentiated sites
tri_plot + annotate('text', x=0.8, y=0.9, label=paste0('Allele frequency threshold: ', alleleFreqThresh))+ 
  annotate('text', x=0.8, y=0.8, label=paste0('Number of sites: ', nrow(GHT_vcf_diff@fix)))

# Repeat the process, but using a higher allele frequency threshold (0.9)
alleleFreqThresh <- 0.9
GHT_vcf_diff <- 
  alleleFreqDiff(vcfR=GHT_vcf, pm=GHT_popmap, p1="GRAV", p2="HYPO", difference=alleleFreqThresh)
GHT_hybridIndex <- hybridIndex(vcfR=GHT_vcf_diff, pm=GHT_popmap, p1="GRAV", p2="HYPO")
# Generate a triangle plot
tri_plot <- triangulaR::triangle.plot(GHT_hybridIndex, colors=GHT_cols, cex=3)
tri_plot + annotate('text', x=0.8, y=0.9, label=paste0('Allele frequency threshold: ', alleleFreqThresh)) + 
  annotate('text', x=0.8, y=0.8, label=paste0('Number of sites: ', nrow(GHT_vcf_diff@fix)))

# Repeat the process, but using a higher allele frequency threshold (1.0)
alleleFreqThresh <- 1
GHT_vcf_diff <- 
  alleleFreqDiff(vcfR=GHT_vcf, pm=GHT_popmap, p1="GRAV", p2="HYPO", difference=alleleFreqThresh)
GHT_hybridIndex <- hybridIndex(vcfR=GHT_vcf_diff, pm=GHT_popmap, p1="GRAV", p2="HYPO")
# Generate a triangle plot
tri_plot <- triangulaR::triangle.plot(GHT_hybridIndex, colors=GHT_cols, cex=3)
tri_plot + annotate('text', x=0.8, y=0.9, label=paste0('Allele frequency threshold: ', alleleFreqThresh))+ 
  annotate('text', x=0.8, y=0.8, label=paste0('Number of sites: ', nrow(GHT_vcf_diff@fix)))

# TRIANGLE PLOTS, FOR CONFERENCE PRESENTATION ----
# Building triangle plots, which are plots of interclass heterozygosity versus hybrid indices,
# and are meant to help distinguish admixture patterns due to hybridization from isolation by distance.
library(triangulaR)
library(vcfR)

# Read in the VCF file generated for 
GHT_vcf <- read.vcfR('populations.snps.vcf')
# Read in a popmap, which is a data.frame with two columns: 'id' and 'pop'
GHT_popmap <- read.table('GravHypoTar_popList.csv', header = TRUE, sep=',')
# Update names in popmap
# colnames(GHT_popmap) <- c('id', 'Taxa')
# GHT_popmap$id <- gsub('GRAV','gravesii',GHT_popmap$id)
# GHT_popmap$id <- gsub('HYPO','hypoleucoides',GHT_popmap$id)
# GHT_popmap$id <- gsub('TARD','tardifolia',GHT_popmap$id)
# Specify a value for the allele frequency difference threshold (only alleles with frequency differences
# of this amount or greater in parental groups are considered for calculating interclass heterozygosity)
alleleFreqThresh <- 0.9
# Create a new vcfR object composed only of sites above the given allele frequency difference threshold
GHT_vcf_diff <- 
  alleleFreqDiff(vcfR=GHT_vcf, pm=GHT_popmap, p1="GRAV", p2="HYPO", difference=alleleFreqThresh)
# GHT_vcf_diff <- 
#   alleleFreqDiff(vcfR=GHT_vcf, pm=GHT_popmap, p1="gravesii", p2="hypoleucoides", difference=alleleFreqThresh)
# Calculate hybrid index and heterozygosity for each sample. Values are returned in a data.frame
GHT_hybridIndex <- hybridIndex(vcfR=GHT_vcf_diff, pm=GHT_popmap, p1="GRAV", p2="HYPO")
# GHT_hybridIndex <- hybridIndex(vcfR=GHT_vcf_diff, pm=GHT_popmap, p1="gravesii", p2="hypoleucoides")

# PLOTTING
# Specify a vector of colors
GHT_cols <- c("#af8dc3", "#7fbf7b", "#bababa", "#878787", "#762a83", "#1b7837")
# Generate a triangle plot
tri_plot <- triangulaR::triangle.plot(GHT_hybridIndex, colors=GHT_cols, cex=3)
# Print the allele frequency threshold value and the number of differentiated sites
tri_plot + annotate('text', x=0.8, y=0.9, label=paste0('Allele frequency threshold: ', alleleFreqThresh)) + 
  annotate('text', x=0.8, y=0.8, label=paste0('Number of sites: ', nrow(GHT_vcf_diff@fix))) +
  guides(color=guide_legend(title="Taxa")) + theme(legend.position = c(0.5,0.3))
