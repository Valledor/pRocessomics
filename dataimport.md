# pRocessomics
pRocessomics makes the integration of different omic levels easy. In this page you will learn how to prepare your data and import it to be used in pRocessomics.

## Data import
The first problem that new users face when using bioinformatic tools is introducing their data in the script/program/... To ease this problem pRocessomics has a built-in function to import data directly from an excel file. This function is called **importfromexcel**. This function creates a dataset list to be further used in pRocessomics, and also perform some initial quality and consistency checks. If data is not OK it will give the user some hints to correct it.

The best/quickest approach to prepare your data to be used in pRocessomics is download a template or an example file and prepare your data in the same way. 
* [Download template file](/datasets/datasettemplate.xlsx)
* [Download UV dataset](/datasets/uvdataset.xlsx)

Please consider that if you want to use this parser the different omic levels within your dataset (from 1 to n) must be integrated in the same Excel file, with one sheet per level. Each Excel sheet should be named with the name you want to give to each omic level. 
As a golden rule: cases (samples) data will be in rows, while variables will be in columns. First column must indicate the name of the cases or samples, these names should be the same across the different omic levels. All omic levels should have the same number of cases (this is a current limitation of this package). We recommend that second and third column will correspond to the different treatments you applied to your data. Most of the times you will only need one column, since only one treatment/factor was tested. Inmmediately after these columns, variables (proteins, metabolites, rnas, ...) should be indicated. We strongly recommend that you employ protein, metabolite, gene, ... accessions instead of defline or description for naming variables. Later you will be able to upload an annotation matrix to get the most of your data.

With an Excel file properly formated, importing data is trivial since we have developed an interactive wizard that does the job. First run the function. It has only one argument, the name you want to give to your dataset, but if it is not present it will ask you for one.

```
  > importfromexcel()
Please provide a name for your dataset (not numeric: 123 -> wrong; 123a -> good) and press [enter]: mydataset
```
After pressing enter a dialog window will appear to choose the file with the dataset. Select the file you want and press "Open". Please **do not** use double click for selecting the file. Just afterwards importfromexcel() will read your file and show you a portion to define some useful information to characterize the structure of your table, and also perform some checks to ensure data structure is OK:
```
Below you will find a small portion of your dataset so you can easily
check which is your starting row, starting column, and treatment
columns.
          Col 1     Col 2             Col 3             Col 4             Col 5
         Sample Treatment           6093823         205829383         383167485
Row 1 Control-1         C 0.330974271602681 0.551627114553486 0.288355529976485
Row 2 Control-2         C 0.298562088023843  0.39808278403179 0.346822251800819
Row 3 Control-3         C  0.24035928562732 0.437014137187948 0.304592017565927
Row 4    1d2h-1        2h 0.310704181881698 0.414261964840912 0.144370402662572
Row 5    1d2h-2        2h 0.292817011254695 0.390419263868094 0.136052556983632
Enter the number of the column containing sample names and press [enter]: 1
Enter the number of the column with the first treatment and press [enter]: 2
Enter the number of the column with the second treatment and press [enter] (optional, press [enter] if you only have one treatment/level): 
Enter the number of the column with the first variable (numerical, not treatment) and press [enter]:  3
Enter the number of the row with the first case (individual) and press [enter]: 1


Removing empty columns.

Job finished!
```
Once tables are read, importfromexcel() will automatically remove empty columns (variable name but no data), and check if rows are also blank. Function will also remove all rows from 1 to the row with the first case. If everything is OK you should get `Job Finished!` message. Furthermore a file with information about experiment will be created with the name of the dataset and the suffix `_dsettings`

If something goes wrong you should will get a detailed error message with some suggested solutions:
```
TABLE ERROR: Error in case names in proteins metabolites dataset(s).

Hint: Case (sample) names should be indicated in the first column of each tables, futhermore these names should be the same across datasets. Cases should also have the same order in all sheets.
```

Obviously you can also do that in the hard way. In this case you can import the different tables of your different datasets creating a list. List should be named. It is not necessary to create `_dsettings`file by yourself, it will be created when using `preprocessing_wizard()`.
```
  > mydatalist <- list("proteins"=table1,"metabolites"=table2,...)
```
