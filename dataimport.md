# pRocessomics
pRocessomics makes the integration of different omic levels easy. In this page you will learn how to prepare your data and import it to be used in pRocessomics.

## Data import
The first problem that new users face when using bioinformatic tools is introducing their data in the script/program/... To ease this problem pRocessomics has a built-in function to import data directly from an excel file. This function is called **importfromexcel**. This function creates a dataset list to be further used in pRocessomics, and also perform some initial quality and consistency checks. If data is not OK it will give the user some hints to correct it.

The best/quickest approach to prepare your data to be used in pRocessomics is download a template or an example file and prepare your data in the same way. 
* [Download template file](/datasets/datasettemplate.xlsx)
* [Download UV dataset](/datasets/uvdataset.xlsx)

Please consider that if you want to use this parser the different omic levels within your dataset (from 1 to n) must be integrated in the same Excel file, with one sheet per level. Each Excel sheet should be named with the name you want to give to each omic level. 
As a golden rule: cases (samples) data will be in rows, while variables will be in columns. First column must indicate the name of the cases or samples, these names should be the same across the different omic levels. All omic levels should have the same number of cases (this is a current limitation of this package). We recommend that second and third column will correspond to the different treatments you applied to your data. Most of the times you will only need one column, since only one treatment/factor was tested. Inmmediately after these columns, variables (proteins, metabolites, rnas, ...) should be indicated. We strongly recommend that you employ protein, metabolite, gene, ... accessions instead of defline or description for naming variables. Later you will be able to upload an annotation matrix to get the most of your data.

With an Excel file properly formated, importing data is trivial:

```
  > mydata <- datalist <- importfromexcel(excelfilename = "uvdataset.xlsx",treatment1col = 2, initialrow = 1,initialcol = 3)
```
The function has different arguments:

  `excelfilename`   xls or xlsx file with all omic levels to be analyzed. Each dataset should be in a different sheet. Dataset will be named after sheet names. If no name is provided, a dialog box for selecting file will appear.

  `treatment1col`   Column number in which treatment 1 is indicated.

  `treatment2col`   (optional) Column number in which treatment 2 is indicated.

  `initialrow`  First row of numeric data within datasets. Non significant rows (from first row to initialrow) will be removed.

  `initialcolumn`   Column of the first variable of the dataset.

Please check package documentation for more/updated information about the different arguments you must use.


Obviously you can also do that in the hard way. In this case you can import the different tables of your different datasets
```
  > mydatalist <- list("proteins"=table1,"metabolites"=table2,...)
```
