# First step
Download and install this package from repository

# Data import
Already attached dataset is a list class object and contains four omic layers: proteome, metabolome, gene expression and, physiology. This dataset comes from an adaptive study of pinus radiata upon UV stress (Pascual et al., 2016).
Please note all listed omic layers are named and follow the same main structure, the treatment column is always located at the same position and also do the numeric data (columns & rows)

**Having an adequate data structure is required for working with pRocessomics**. If you overcome this obstacle, the use of this package will be very simple.

This package is meant to work with lists, even if you only want to process a single omic dataset you need your table to be in a list.
The different elements in the list should be named (proteins, metabolites, etc.) and have the same number of rows. Each row correspond to a sample. Samples should have the same order within the different datasets to be analyzed (if you want to integrate the different datasets)
If you have troubles importing your data into pRocessomics, **please check this tutorial**

For this case stady you only need to load sample dataset 
