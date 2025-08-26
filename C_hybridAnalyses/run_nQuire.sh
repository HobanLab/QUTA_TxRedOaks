#!/bin/bash

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% ASSESSING POLYPLOIDY USING NQUIRE %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# This script calls the nQuire commands used to assess ploidy levels in the QUTA/Texas red oak samples.
# For each .bam file in the input directory, the following steps are taken:
#
# 1. Create the .bin file, using the -x flag to store reference IDs and positional information (genomic coordinates)
# 2. Denoise the base frequency histogram, scaling it down to more easily detect peaks in the histogram of base frequencies
#
# Note: the coverage threshold (-c argument) is externally passed as the first argument to the script
#
# Once denoised .bin files are created, the lrdmodel command is used to describe the histogram of each file
# as a mixture of Gaussian distributions with varying means and mixture proportions. The Expectation-Maximization (EM)
# algorithm is then run, and this is used to calculate the delta log-likelihood for each model (diploid, triploid, tetraploid).
# The lowest delta log-likelihood value indicates the best supported ploidy scenario.

# The output from lrdmodel contains 8 tab-separated columns:
# 1. Filename
# 2. Free model maximized log-likelihood
# 3. Diploid fixed model maximized log-likelihood
# 4. Triploid fixed model maximized log-likelihood
# 5. Tetraploid fixed model maximized log-likelihood
# 6. Diploid delta log-likelihood
# 7. Triploid delta log-likelihood
# 8. Tetraploid delta log-likelihood
#
# !!! NOTE: THIS SCRIPT NEEDS TO BE RUN FROM THE DIRECTORY IN WHICH THE .BIN FILES ARE CREATED !!!
#
# %%% VARIABLES
# Specify a filepath variable for the directory where input filtered sequences are stored
inputDir=/RAID1/QUTA_TX_RedOaks/Genotyping/Alignment/BWA2/SE_Only
# Declare a variable capturing all of the .bam files in the directory; these will be looped over
inputFiles="$inputDir/*.bam"
# Specify number of threads for building models
n_threads=16
# Declare the coverage argument passed by the user to the script
echo "%%% Minimum coverage argument: $1x %%%"

# Loop through the samples listed in the input files variable
for samplePath in $inputFiles
do
        # Create a variable containing strictly the sample name (without the leading filepath or file extension)
        sample="${samplePath##*/}"
	sample=${sample%.*}
        # Create sample name variables for input and output files
        bam=${inputDir}/${sample}.bam
        # Create the binary file, using the -x flag to store reference IDs and positional information (genomic coordinates)
	# The filtering defaults (-f: minimum frequency 0.2; -q: minimum mapping quality 1) are used
	# Coverage threshold (-c) is user-defined; if none is specified, default value (10x) is used
	nQuire create -b $bam -o ${sample} -x -c $1
        # Denoise the binary file, to moreeasily detect peaks in the histogram of base frequencies
	nQuire denoise ${sample}.bin -o ${sample}_d
done

# Build the GMMs for each denoised .bin file. Default parameters are used.
nQuire lrdmodel -t $n_threads *_d.bin
