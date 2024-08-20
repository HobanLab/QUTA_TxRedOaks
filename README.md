# Overview
This repository contains the code used for the population genetic analyses of _Quercus tardifolia_ and other Texas red oaks. The goals of these analyses are:
1. determining how genetically distinctive _Q. tardifolia_ (QUTA) is
2. understanding the taxonomic entities which exist in this group
3. assess any admixture patterns in other taxonomic names in this group (_robusta_, _graciliformis_, etc.)

The genetic data for these analyses came from the [Hipp Lab](http://systematics.mortonarb.org/lab/) at the Morton Arboretum. The code in this repository is exploratory and in development. 

## Repository layout
There are three primary folders in this repository.

1. [`Genotyping`](https://github.com/HobanLab/QUTA_TRO/tree/main/Genotyping) contains the BASH script and input files used to perform demultiplexing (when necessary), QC, reference alignment, and genotyping. There's also 1 R script, which was used
to assess the level of polymorphic error between replicate samples.

2. [`Clustering`](https://github.com/HobanLab/QUTA_TRO/tree/main/Clustering) contains the scripts used to run STRUCTURE, build neighbor-joining trees, and perform discriminant analysis of principal components. This folder is divided into subfolders
based on which clustering step (prelim, Clust2, or Clust3) the script corresponds to. The outputs from the DAPC runs are also stored in this folder.

3. [`hybridAnalyses`](https://github.com/HobanLab/QUTA_TRO/tree/main/hybridAnalyses) contains the scripts used to validate the hybridization patterns found during the clustering steps. As of right now, this strictly consists of calculation of f3 statistics--
but this will likely change over time.
