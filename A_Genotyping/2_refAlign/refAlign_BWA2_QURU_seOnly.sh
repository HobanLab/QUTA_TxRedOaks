#!/bin/bash

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%% CREATE ALIGNMENTS WITH QUERCUS RUBRA REFRENCE GENOME %%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# This script aligns demultiplexed and filtered paired-end fastq files (processed using the Stacks process_radtags command) to the Quercus rubra reference genome
# (Kapoor et al. 2023) using the BWA software. It does so in a loop, which loops over a barcode file that lists the names of each sample processed.
# Additionally, it uses the samtools flagstat command to generate alignment statistics for each sample, and write these to a summary file. After paired-end reads,
# single-end reads are aligned. The samples processed in this script are a subset of those processed in the demultiplexing/QC steps, as barcode files only select
# the samples that are Texas (and in a few cases, Mexican) red oaks.

# Filepath variables are specified to direct to the appropriate folders on the RAID1 drive. This script needs to be run from the directory where the Q. rubra reference
# genome is stored (in order for the bwa mem command to not error). The contrast with this script and BWA1 is that the -P flag is not included, for paired-end reads
# (as this was found to cause the properly paired percentage of reads to be zero). In this script, for both paired-end and single-end samples, only forward reads are
# processed (in order to make sure all samples can be included in the same library).

# %%% PAIRED-END READS %%%%
echo "%%% PROCESSING PAIRED-END SAMPLES %%%"

# %%% FILEPATH VARIABLES
# Specify a filepath variable for the directory where input filtered sequences are stored
inputDir=/RAID1/QUTA_TX_RedOaks/Genotyping/Demux_QC/Output
# Specify a filepath variable for the directory where barcode files are stored. Note these are different than the barcode files used for demultiplexing/QC,
# because only red oak samples are included (rather than all samples included on a plate)
barcodeDir=/RAID1/QUTA_TX_RedOaks/Genotyping/Alignment/BWA2/Barcodes
# Specify a filepath variable for the directory where results (alignment files) will be stored
outputDir=/RAID1/QUTA_TX_RedOaks/Genotyping/Alignment/BWA2/bam_files

# %%% OUTER LOOP
# Loop over the names of each sample plate, storing the name. This is used to direct to the correct input sequence
for plate in C1322 C1343 C1440
do
	# Print the name of the plate being processed, so users can keep track in the printed outputs
	progMsg="*** PLATE CURRENTLY BEING PROCESSED: $plate ***"
	echo "${progMsg}"
	# %%% INNER LOOP
	# Loop over the sample names listed in the barcode file used for alignment (selecting only the 2nd column, which has the sample names)
	cut -f2 ${barcodeDir}/${plate}_Barcodes.txt |
	while read sample; do
        	# Create filename variables for fastq Read 1 (forward read). Reverse reads are not processed
	        fq1=${inputDir}/${plate}/${sample}.1.fq.gz # Forward read
        	# Create sample name variables for output alignment files
        	sam=${outputDir}/${sample}.sam
        	bam=${outputDir}/${sample}.bam
		# Align reads using BWA mem, referring to the indexed Q. rubra reference genome (using database name specified in bwa index command)
		# Do not use the -P flag (description: perform SW to rescue missing hits only but do not try to find hits that fit a proper pair.)
	        bwa mem Qrubra_db $fq1 -t 30 -o $sam
	        # Compress alignments, then pass the output onto sort to generate sorted BAM file
        	samtools view -bh --threads 30 $sam | samtools sort --threads 30 -o $bam
	        # Generate basic alignment statistics using flagstat, and output them to a summary file (one file used for all samples)
		echo ${sample} >> ${outputDir}/PE_algnSummary_BWA2.tsv
		samtools flagstat -O tsv $bam >> ${outputDir}/PE_algnSummary_BWA2.tsv
	done
done

# %%% SINGLE-END READS %%%
echo "%%% PROCESSING SINGLE-END SAMPLES %%%"

# %%% FILEPATH VARIABLES
# Specify a filepath variable for the directory where input filtered sequences are stored
inputDir=/RAID1/QUTA_TX_RedOaks/Genotyping/Demux_QC/Output/SingleEndPlate
# Declare a variable capturing all of the .fq.gz files in the paired-end plate directory; these will be looped over
inputFiles="$inputDir/*.fq.gz"
# Specify a filepath variable for the directory where results (alignment files) will be stored. We specify a subdirectory for single-end reads
outputDir=/RAID1/QUTA_TX_RedOaks/Genotyping/Alignment/BWA2/bam_files/SE_Samples

# Loop through the samples listed in the input files variable
for samplePath in $inputFiles
do
        # Create a variable containing strictly the sample name (without the leading filepath)
	sample="${samplePath##*/}"
	# Create sample name variables for output alignment files
        sam=${outputDir}/${sample}.sam
        bam=${outputDir}/${sample}.bam
        # Align reads using BWA mem, referring to the indexed Q. rubra reference genome (using database name specified in bwa index command)
	bwa mem -P Qrubra_db $samplePath -t 30 -o $sam
        # Compress alignments, then pass the output onto sort to generate sorted BAM file
        samtools view -bh --threads 30 $sam | samtools sort --threads 30 -o $bam
        # Generate basic alignment statistics using flagstat, and output them to a summary file
        echo ${sample} >> ${outputDir}/SE_algnSummary_BWA2.tsv
        samtools flagstat -O tsv $bam >> ${outputDir}/SE_algnSummary_BWA2.tsv
done
