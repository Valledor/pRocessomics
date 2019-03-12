# pRocessomics
pRocessomics makes the integration of different omic levels easy. How to import your data.

## Data import
The first problem that new users face when using bioinformatic tools is introducing their data in the script/program/... To ease this problem pRocessomics has a built-in function to import data directly from an excel file. This function is called **importfromexcel**. This function creates a dataset list to be further used in pRocessomics, and also perform some initial quality and consistency checks. If data is not OK it will give the user some hints to correct it.

The best/quickest approach to prepare your data to be used in pRocessomics is download a template or an example file and prepare your data in the same way. 

Please consider that if you want to use this parser the different omic levels within your dataset (from 1 to n) must be integrated in the same Excel file, with one sheet per level. Each Excel sheet should be named with the name you want to give to each omic level. 
As a golden rule: cases (samples) data will be in rows, while variables will be in columns. First column must indicate the name of the cases or samples, these names should be the same across the different omic levels. All omic levels should have the same number of cases (this is a current limitation of this package). We recommend that second and third column will correspond to the different treatments you applied to your data. Most of the times you will only need one column, since only one treatment/factor was tested. Inmmediately after these columns, variables (proteins, metabolites, rnas, ...) should be indicated. We strongly recommend that you employ protein, metabolite, gene, ... accessions instead of defline or description for naming variables. Later you can upload an annotation matrix to get the most of your data.



```
  > mydata <- importfromexceldatalist(excelfilename=NULL,treatment1col=NULL,treatment2col=NULL, initialrow=NULL, initialcol=NULL)
```

Obviously you can also do that in the hard way. In this case you can import the different tables of your different datasets
```
  > mydatalist <- list("proteins"=table1,"metabolites"=table2)
```
