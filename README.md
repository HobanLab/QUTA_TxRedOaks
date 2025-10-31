# Overview
This repository contains the code used for the population genetic analyses of _Quercus tardifolia_ and other Texas red oaks. The goals of these analyses are:
1. determining how genetically distinctive _Q. tardifolia_ (QUTA) is
2. understanding the taxonomic entities which exist in this group
3. assess any admixture patterns in other taxonomic names in this group (_robusta_, _graciliformis_, etc.)

The genetic data for these analyses came from the [Hipp Lab](http://systematics.mortonarb.org/lab/) at the [Morton Arboretum](https://mortonarb.org/science/center-for-tree-science/). The code in this repository is exploratory and in development. 

## Repository layout
There are four primary folders in this repository (ordered alphabetically, in approximate processing order):

1. [`Genotyping`](https://github.com/HobanLab/QUTA_TRO/tree/main/A_Genotyping) contains the BASH script and input files used to perform demultiplexing (when necessary), QC, reference alignment, and genotyping. There's also 1 R script, which was used
to assess the level of polymorphic error between replicate samples.

2. [`Clustering`](https://github.com/HobanLab/QUTA_TRO/tree/main/B_Clustering) contains the scripts used to run [STRUCTURE](https://web.stanford.edu/group/pritchardlab/structure.html), build neighbor-joining trees, and perform discriminant analysis of principal components (DAPC) and principal components analysis (PCA). This folder is divided into subfolders
based on which clustering step (prelim, Clust2, or Clust3), as well as a folder for DAPC outputs. Note that the R scripts used for visualizing the STRUCTURE results are contained in this folder.

3. [`hybridAnalyses`](https://github.com/HobanLab/QUTA_TRO/tree/main/C_hybridAnalyses) contains the scripts used to validate the hybridization patterns found during the clustering steps. There are two subfolders: one containing the scripts used to run the species networks applying quartets ([SNaQ](https://juliaphylo.github.io/SNaQ.jl/stable/)) software, and one containing the R script and output of the f-stats analysis (which presently isn't included in the main results of this project.) 

  There is also an R script used to run the [nQuire](https://github.com/clwgg/nQuire) software (for   assessing polyploidy) and for conducting heterozygosity analyses (e.g. using   [triangulaR](https://omys-omics.github.io/triangulaR/index.html)).

4. [`Mapping`](https://github.com/HobanLab/QUTA_TxRedOaks/tree/main/D_Mapping) contains the R script used to generate an HTML map of the samples. This map was largely used internally.
