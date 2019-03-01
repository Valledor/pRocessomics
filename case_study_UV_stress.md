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

# Stage 1, data pre-processing
Before exploring our data, few considerations may be taken into account: as omics datasets are often obtained from mass spectrometry techniques, such as proteome and metabolome, some values can be missed, or may need a global abundance scaling. To overcome these issues, the first step in our analysis is to run preprocess_omic_list function. This function will perform missing data imputation, empty columns removal and abundance balancing.

## Missing value imputation
There are three available options: using a RandomForest algorithm, k-nearest neighbour and do not imput missing values. This argument can be provided as a vector defining the desired operation for each omic layer in the dataset. In our case study we will employ different methods to illustrate the capability of pRocessomics

## Abundance balancing
pRocessomics has implemented three kinds of abundance scaling: sample (relative amount of each variable within samples, as percentage; gives really small values, not easy to interpret), treatment intensity (normalized by all samples within each treatment) and average intensity (percentage multiplied for the average total intensity of samples, this is proportional to the first approach, but numbers are bigger, in "real-world" ranges, and in consequence easier to interpret). When dealing with biological samples it is very important to make sure that numerical values can be compared across samples. We must control differences due to bad sample loading, missquantification of initial amounts, etc. 

### Source code
```
  > my.preprocessed.list <- preprocess_omic_list(datalist = datalist, initialrow = 1,
      initialcolumn = 2,treatment1col = 1,treatment2col = NULL,treatment = 1, 
      imputation = c("KNN","RF","none","none"),abdbal = c("AvgIntensity","AvgIntensity","none","none"),
      threshold = 0.34,k = 2,parallel = TRUE)
```
  Expected output:
```
  MISSING VALUE IMPUTATION
Multiple processor cores (4) will be used. It may take a while...
KNN method will be used for proteome dataset
RF method will be used for metabolome dataset
none method will be used for genexpression dataset
none method will be used for physiology dataset

Processing proteome dataset
Original data contained 5359 zeroes and 0 NAs in  669  and  0  variables, respectively. After processing 796 values present in  439 variables have been considered suitable for imputation according to the defined 0.34 threshold.
Missing values have been imputed according to KNN method.

Processing metabolome dataset
Original data contained 8 zeroes and 0 NAs in  5  and  0  variables, respectively. After processing 5 values present in  5 variables have been considered suitable for imputation according to the defined 0.34 threshold.
  missForest iteration 1 in progress...done!
  missForest iteration 2 in progress...done!
  missForest iteration 3 in progress...done!
  missForest iteration 4 in progress...done!
Missing values have been imputed according to RF method.

Processing genexpression dataset
Original data contained 0 zeroes and 0 NAs in  0  and  0  variables, respectively. After processing no variable has been considered suitable for imputation according to the defined 0.34 threshold.

Processing physiology dataset
Original data contained 0 zeroes and 0 NAs in  0  and  0  variables, respectively. After processing no variable has been considered suitable for imputation according to the defined 0.34 threshold.

REMOVING EMPTY COLUMNS OF ALL DATASETS
Single processor core will be used. It may take a while...

ABUNDANCE BALANCING
Single processor core will be used. It may take a while...
AvgIntensity balancing will be used for proteome dataset
AvgIntensity balancing will be used for metabolome dataset
none balancing will be used for genexpression dataset
none balancing will be used for physiology dataset

Processing proteome dataset
Processing metabolome dataset
Processing genexpression dataset
Processing physiology dataset


SUMMARY: Dataset imputation and balancing
-----------------------------------------
proteome, metabolome, genexpression, physiology datasets were considered
 
 
Missing value imputation.........
proteome, metabolome were respectively imputed employing KNN, RF methods. 
2 neighbors were employed for KNN calculations 
Minimum percentage of significant values per variable/treatment to allow 
imputation was  0.34 
Run in parallel TRUE 
 
Abundance balancing.........
proteome, metabolome were respectively balanced employing AvgIntensity, 
AvgIntensity methods.

 Job finished!
```
## Data transformation and filtering

Sometimes, it may be useful apply a mathematical transformation to the data, in order to get a normal distribution, or eliminate noisy or almost empty variables, to this end transformandselect function can be used as follows:

### Source code
```
> my.transformandselect.list<-transformandselect(datalist = my.preprocessed.list,initialrow = 1,initialcolumn = 2,treatment1col = 1,treatment2col = 1,treatment = 1,transf = c("Log10","Log10","none","z"),varsel = TRUE,varselthld = 0.4,varcoef = FALSE,varcoefthld = NULL)
```
  Expected output:
```
TRANSFORMATION OF DATASETS
Single processor core will be used. It may take a while...
Log10 transformation method will be used for proteome dataset
Log10 transformation method will be used for metabolome dataset
none transformation method will be used for genexpression dataset
z transformation method will be used for physiology dataset

Processing proteome dataset
Processing metabolome dataset
Processing genexpression dataset
Processing physiology dataset
SELECTING VARIABLES BASED ON CONSISTENCY
Single processor core will be used. It may take a while...

Processing proteome dataset
Variable Selection based on consistency. Variable was present at least on 6  cases, or all replicates of a treatment
Initial Variables:  1228  Selected Variables:  1228  Removed Variables:  0
Processing metabolome dataset
Variable Selection based on consistency. Variable was present at least on 6  cases, or all replicates of a treatment
Initial Variables:  118  Selected Variables:  118  Removed Variables:  0
Processing genexpression dataset
Variable Selection based on consistency. Variable was present at least on 6  cases, or all replicates of a treatment
Initial Variables:  6  Selected Variables:  6  Removed Variables:  0
Processing physiology dataset
Variable Selection based on consistency. Variable was present at least on 6  cases, or all replicates of a treatment
Initial Variables:  5  Selected Variables:  5  Removed Variables:  0
done!
```


# Stage 2, Univariate Analysis
Now, we have imputed, balance and, transform our data, we can explore our data. Following a classic approach we will perform an ANOVA test followed by Tukey HSD post hoc. We will employ irradiation time (column 1) as variable to define the different treatments. q-values will be also stimated
```
> Univariate.list<-univariate(datalist = my.transformandselect.list,initialrow = 1,initialcolumn = 2,treatment1col = 1,treatment2col = NULL,treatment = 1,parametric = TRUE,posthoc = TRUE,FDR = TRUE,round = 5)
```
The output of this test will be a list, containing original values, average per treatment, SD per treatment, p and q values. Results can be easily printed on screen or exported into Excel file. I.e. ANOVA p and q values, and also Tukey's for all pairs comparisons for the six first proteins in our dataset:
```
> print(Univariate.list$pvalue$proteome[1:6,1:6])
        X6093823 X205829383 X383167485 X356997196 X383149096 X1168580
p-value  0.01640    0.07305      0e+00    0.00001    0.00000  0.00000
q-value  0.02111    0.08175      0e+00    0.00002    0.00000  0.00000
2h-2d    0.92493    0.96300      0e+00    0.29221    0.00000  0.30820
8h-2d    0.99862    0.99998      7e-05    0.00059    0.99919  0.00045
C-2d     0.97202    0.93936      9e-05    0.99995    0.00000  1.00000
r-2d     0.04394    0.16696      2e-05    0.00635    0.00000  0.00084
```
Table with mean + SD values is also generated
```
> print(Univariate.list$meansd$proteome[1:6,1:6])
           Mean-2d        SD-2d     Mean-2h        SD-2h    
X6093823   "-0.50702" "±" "0.06887" "-0.45496" "±" "0.1092" 
X205829383 "-0.36619" "±" "0.02514" "-0.38547" "±" "0.01758"
X383167485 "0"        "±" "0"       "-0.84328" "±" "0.01759"
X356997196 "-0.02642" "±" "0.01776" "-0.07838" "±" "0.04301"
X383149096 "-0.94231" "±" "0.01936" "0"        "±" "0"      
X1168580   "-0.06125" "±" "0.01942" "-0.10329" "±" "0.04661"
```
Export table with all results of the analysis
```
export_table(Univariate.list, "univariate.xlsx")
