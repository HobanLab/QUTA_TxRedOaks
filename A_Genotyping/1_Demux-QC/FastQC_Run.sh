#!/bin/bash

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% FastQC Script: Demultiplexed Data %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Each line below calls FastQC for a different directory of the demultiplexed data (corresponding to different plates)
# zcat is used to extract the contents of the gzipped fastq files while keeping the files as is
# -o specifies the output directory; -t specifies the number of threads to use

# Paired-end reads
zcat C1322/*.fq.gz | fastqc stdin:C1322 -t 8
zcat C1343/*.fq.gz | fastqc stdin:C1343 -t 8
zcat C1440/*.fq.gz | fastqc stdin:C1440 -t 8
# Single-end reads
zcat SingleEndPlate/*.fq.gz | fastqc stdin:SingleEndPlate -t 8
