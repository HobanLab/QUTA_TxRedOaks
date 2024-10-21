# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% QUTA & TEXAS RED OAKS : SNPS2CF SCRIPT: SPECIES TIPS %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# This is script used to convert a phylip file (or potentially a VCF) to a table of 
# concordance factors (CF), which can then be used by PhyloTools (SNaQ) to generate a phylogenetic
# network. The CF table is saved as a CSV in the end of the script: it contains concordance factor
# values for every possible quartet of individuals included in the input VCF file.

# This script is used to build the CF table of "all possible taxa": the table which includes the individuals in 
# the "base phylogenetic network", but then also all "tardifolia" individuals. The only difference in
# this script versus the AllTips script is the Imap file: for this run, the Imap file groups samples
# into species (rather than keeping them as individuals)

library(pegas)
library(foreach)
library(doMC)
library(parallel)

# %%% VARIABLES %%% ----
# Read in the functions used to convert a phylip file to a CF table
QUTA_SNPs2CF_wd <- '/RAID1/QUTA_TX_RedOaks/HybridAnalyses/SNPs2CF/GH_Repo_SNPs2CF/functions_v1.6.R'
source(QUTA_SNPs2CF_wd)
# Set working directory to folder containing phylip file
QUTA_SNaQ_wd <- '/RAID1/QUTA_TX_RedOaks/HybridAnalyses/SNPs2CF/SpeciesTips/'
setwd(QUTA_SNaQ_wd)
# Set up relevant cores 
num_cores <- detectCores() - 4 

# Specify filepaths to input phylip and map files, and to output CSV
INPUT_VCF <- 'Input/QUTA_Clust3_SNPs2CF_SpeciesTips.vcf'
OUTPUT_phylip <- 'Input/QUTA_Clust3_SNPs2CF_SpeciesTips.phylip'
INPUT_map <- 'Input/QUTA_Clust3_SNPs2CF_SpeciesTips_Imap.txt'
OUTPUT_CFtable <- 'Output/Clust3_SNPs2CF_SpeciesTips_CFtable2.csv'

# %%% CONVERT VCF TO PHYLIP ----
# This command only needs to be run once; once run, a phylip file will be written to disk. This
# phylip file will be the input for the SNPs2CF command (see below). 
# vcf2phylip(vcf.name=INPUT_VCF, total.SNPs=86142, random.phase=T, 
#            output.name=OUTPUT_phylip,  cores=num_cores)
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
