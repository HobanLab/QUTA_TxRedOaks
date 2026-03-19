[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18306919.svg)](https://doi.org/10.5281/zenodo.18306919)

# Overview
This repository contains the code used for the population genetic analyses of _Quercus tardifolia_ and other Texas red oaks. The goals of these analyses are:
1. determining how genetically distinctive _Q. tardifolia_ (QUTA) is
2. ff QUTA is found to be a hybrid, determine the hybrid class of QUTA individuals
3. understanding the taxonomic entities which exist in this group

The genetic data for these analyses came from the [Hipp Lab](http://systematics.mortonarb.org/lab/) at the [Morton Arboretum](https://mortonarb.org/science/center-for-tree-science/). The code in this repository is exploratory and in development. 

## Repository layout
There are four primary folders in this repository (ordered alphabetically, in approximate processing order):

1. [`Genotyping`](https://github.com/HobanLab/QUTA_TRO/tree/main/A_Genotyping) contains the BASH script and input files used to perform demultiplexing (when necessary), QC, reference alignment, and genotyping. There's also one R script, which was used
to assess the level of polymorphic error between replicate samples. Barcode files used during multiplexing are included in the 1_Demux-QC folder.

2. [`Clustering`](https://github.com/HobanLab/QUTA_TRO/tree/main/B_Clustering) contains the scripts used to run [STRUCTURE](https://web.stanford.edu/group/pritchardlab/structure.html), build neighbor-joining trees, and perform discriminant analysis of principal components (DAPC) and principal components analysis (PCA). This folder is divided into subfolders
based on STRUCTURE clustering runs (prelim, Clust2, or Clust3; Clust3 includes the most relevant results to the manuscript), as well as a folder for DAPC outputs. Note that the R scripts used for visualizing the STRUCTURE results are contained in this folder.

3. [`hybridAnalyses`](https://github.com/HobanLab/QUTA_TRO/tree/main/C_hybridAnalyses) contains the scripts used to validate the hybridization patterns found during the clustering steps. There are two subfolders: one containing the scripts used to run the species networks applying quartets ([SNaQ](https://juliaphylo.github.io/SNaQ.jl/stable/)) software, and one containing the R script and output of the F-stats analysis (which presently isn't included in the main results for this project.) There is also an R script used to run the [nQuire](https://github.com/clwgg/nQuire) software (for assessing polyploidy) and for conducting heterozygosity analyses (e.g. using [triangulaR](https://omys-omics.github.io/triangulaR/index.html)). Note that the `GravHypoTar_HzAnalysis.R` script includes the code for the triangulaR plots.

4. [`Mapping`](https://github.com/HobanLab/QUTA_TxRedOaks/tree/main/D_Mapping) contains two R scripts. 
	1. `QUTA_sampleMapping.R` was used to generate an HTML map of all the samples; this map was largely used internally. 
	2. The `GravHypoTar_RangeMapping.R` script was used to build a range map of *Q. gravesii* and *Q. hypoleucoides*, and show sampled individuals on that map.
