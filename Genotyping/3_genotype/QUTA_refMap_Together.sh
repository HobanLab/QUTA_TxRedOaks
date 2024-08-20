#!/bin/bash

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%% CALL SNPS USING STACKS REF_MAP.PL SCRIPT: PROCESS PAIRED-END AND SINGLE-END RESULTS TOGETHER %%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# This script calls the ref_map.pl Stacks command, in order to call scripts from the BWA1 alignment of QUTA/Texas red oak samples with the Quercus rubra reference genome.
# Calling this script runs the gstacks and populations Stacks modules. The gstacks module simply genotype individuals and and phases SNPs at each locus, generating a set
# of haplotypes. The populations module is used to filter loci based on different preferences, and the generate output file formats to pass onto downstream analyses.

# One command is called, which directs to a folder that contains a combination of paired-end and single-end reads, or a folder which contains strictly single-end reads.

# %%% REF_MAP.PL
ref_map.pl -T 28 --samples /RAID1/QUTA_TX_RedOaks/Genotyping/Alignment/BWA2/SE_Only --popmap ./Popmaps/QUTA_All_popmap -o ./Output/All_seReadsOnly
