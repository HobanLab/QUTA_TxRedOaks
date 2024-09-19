# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% QUTA & TEXAS RED OAKS : SNPS2CF SCRIPT %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# This is script used to convert a phylip file (or potentially a VCF) to a table of 
# concordance factors (CF), which can then be used by PhyloTools to generate a phylogenetic
# network.

library(pegas)
library(foreach)
library(doMC)
library(parallel)

# %%% VARIABLES %%% ----
# Read in the functions used to convert a phylip file to a CF table
QUTA_SNPs2CF_wd <- '/RAID1/QUTA_TX_RedOaks/HybridAnalyses/GH_Repo_SNPs2CF/functions_v1.6.R'
source(QUTA_SNPs2CF_wd)
# Set working directory to folder containing phylip file
QUTA_SNaQ_wd <- '/RAID1/QUTA_TX_RedOaks/HybridAnalyses/SNPs2CF/'
setwd(QUTA_SNaQ_wd)
# Set up relevant cores 
num_cores <- detectCores() - 4 

# Specify filepaths to input phylip and map files, and to output CSV
INPUT_VCF <- 'input/QUTA_Clust3_SNaQ1.vcf'
OUTPUT_phylip <- 'input/QUTA_Clust3_SNaQ1.phylip'
INPUT_map <- 'input/QUTA_Clust3_SNaQ1_Imap.txt'
OUTPUT_CFtable <- 'output/Clust3_SNaQ1_CFtable_8sp.csv'

# %%% CONVERT VCF TO PHYLIP ----
# This command only needs to be run once; once run, a phylip file will be written to disk. This
# phylip file will be the input for the SNPs2CF command (see below). 
vcf2phylip(vcf.name=INPUT_VCF, total.SNPs=71483, random.phase=T, 
           output.name=OUTPUT_phylip,  cores=num_cores)
# NOTE: converting the VCF to a phylip file leads to 2 lines per individual in the PHYLIP 
# file to be written (one for each haplotype, labeled _0 and _1). The map file therefore also
# needs to have these two lines per individual specified.

# %%% RUN SNPS2CF COMMAND ----
# Print starting time
startTime <- Sys.time() 
print(paste0('%%% SNPs2CF START: ', startTime))

# In this command, 2 individuals (4 haplotypes) per species are specified in the phylip file. 
# Because the between.sp.only command is set to TRUE, all of the sampled quartets will only 
# involve different species (and no more than 1 individual per species). Bootstrap value
# set to TRUE (default)
output <- SNPs2CF(seqMatrix=OUTPUT_phylip, ImapName=INPUT_map, between.sp.only=TRUE, 
                  outputName=OUTPUT_CFtable, save.progress=FALSE, cores=num_cores)
print(paste0('\n', '%%% Completed running the SNPs2CF command!!! %%%'))

# Calculate and print runtime
endTime <- Sys.time() 
print(paste0('%%% SNPs2CF END: ', endTime))
cat(paste0('\n', '%%% TOTAL RUNTIME: ', endTime-startTime))

# # TEST COMMAND: commented out, but used as an example
# output <- SNPs2CF(seqMatrix=OUTPUT_phylip, ImapName=INPUT_map, between.sp.only = TRUE, 
#                   max.SNPs = 50, bootstrap=FALSE, outputName=OUTPUT_CFtable, save.progress=FALSE)

# # %%% ARCHIVE: ERRORS WHEN READING IN PHYLIP FILES ----
# # There were two kinds of phylip files generated using the Stacks populations module:
# # 1. SNaQ1: specifies all 8 samples belonging to the same population (termed 'wild')
# # 2. SNaQ2: specifies all 8 samples belonging to the unique populations (shortened sample name)
# # Option #1 generated a phylip files where all samples have the name 'wild', which was not 
# # usable by SNPs2CF. Therefore, SNaQ2 is analyzed here.
# 
# # %%% RUN SNPS2CF COMMAND: SNaQ2
# INPUT_phylip <- 'input/Stacks_phylips/QUTA_Clust2_SNaQ3.phylip'
# INPUT_map <- 'input/Stacks_phylips/QUTA_Clust3_SNaQ2_Imap.txt'
# OUTPUT_CFtable <- 'output/Clust3_SNaQ2_CFtable.csv'
# 
# # Run SNPs2CF command
# output <- SNPs2CF(seqMatrix=INPUT_phylip, ImapName=INPUT_map, between.sp.only = TRUE, 
#                   max.SNPs = 25, bootstrap=FALSE, outputName=OUTPUT_CFtable, save.progress=FALSE)
# # Error: "Error in SNP matrix. Each individual has to have a unique name."
# 
# # Removing spaces after some of the sample names in the phylip file (361 and 410; SNaQ3) does 
# # not solve the issue.
