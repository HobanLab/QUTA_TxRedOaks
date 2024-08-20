#!/bin/bash

# Single line command which creates a "whitelist" text file, which specifies a random collection of 10,000 loci to pass along to the Stacks populations module
grep -v "^#" pop_R98/populations.sumstats.tsv | cut -f 1 | sort | uniq | shuf | head -n 10000 | sort -n > QUTA_R98_WL_10k.txt

# This command comes from the Stacks manual, found here: https://catchenlab.life.illinois.edu/stacks/comp/populations.php
# A description of the steps involved in this command are provided below

# 1. grep pulls out all the lines in the sumstats file, minus the commented header lines. The sumstats file contains all the polymorphic loci in the analysis.
# 2. cut out the second column, which contains locus IDs
# 3. sort those IDs
# 4. reduce them to a unique list of IDs (remove duplicate entries)
# 5. randomly shuffle those lines
# 6. take the first 10,000 of the randomly shuffled lines
# 7. sort them again and capture them into a file

# The populations.sumstats.tsv file referenced should direct to the original populations dataset.
