#!/bin/bash/

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% QUTA TEXAS RED OAKS: RUN STRUCTURE %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# This script initiates the STRUCTURE runs used to cluster individuals. It refers to a .str file generated using the Stacks populations module.
# Two loops are used to run STRUCTURE. The outer loop iterates through the K values explored: here, 2 to 12
# The inner loop is used to run STRUCTURE in replicate 20 times, for each K value.
# The mainparams.txt file referred to in the STRUCTURE call specifies the number of loci and individuals to analyze.

# To get the start processing time, get the current system time, taken at the beginning of the run, and print it
now=$(date)
echo "Current date/time: $now"
# End processing time can be taken from the output file timestamps (i.e. nohup.out)

# Processing the STRUCTURE file: Stacks prints STRUCTURE file with a populations column, which causes issues when pusing the file through STRUCTRE
# To remove this column (the second column), use the cut command and create a new file
# -f2 specifies cutting the 2nd column, and complement implies keep everything EXCEPT that 2nd column
cut -f2 --complement QUTA_TRO_Clust2_R98_WL10k.str > QUTA_TRO_Clust2_R98_WL10k_F.str

# STRUCTURE loop
# For K values from 2 - 7...
for k in {2..7} ;
do
    # run STRUCTURE 20 times (20 replicates per K value)
    for r in {1..20} ;
    do
	$STR_PATH -i QUTA_TRO_Clust2_R98_WL10k_F.str -m mainparams.txt -e extraparams.txt -K $k -o Output/QUTA_TRO_Clust2_$k-$r > Output/Outfiles/outfile_$k_$r.txt  &
	# Sleep command is to prevent sending all STRUCTURE calls simultaneously
        sleep 3s
    done
done
echo "All runs started!"
