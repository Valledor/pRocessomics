#case study 1
install.packages("pRocessomics")
library("pRocessomics")

data(datalist)

names(datalist)

print(datalist$proteome[1:10,1:4])
print(datalist$metabolome[1:10,1:4])
print(datalist$genexpression[1:10,1:4])
print(datalist$physiology[1:10,1:4])

my.preprocessed.list<-preprocess_omic_list(datalist = datalist,initialrow = 1,initialcolumn = 2,treatment1col = 1,treatment2col = NULL,treatment = 1,imputation = c("KNN","RF","none","none"),abdbal = c("AvgIntensity","AvgIntensity","none","none"),threshold = 0.34,k = 2,parallel = TRUE )
#MISSING VALUE IMPUTATION
#Multiple processor cores (4) will be used. It may take a while...
#KNN method will be used for proteome dataset
#RF method will be used for metabolome dataset
#none method will be used for genexpression dataset
#none method will be used for physiology dataset

#Processing proteome dataset
#Original data contained 5359 zeroes and 0 NAs in  669  and  0  variables, respectively. After processing 796 values present in  439 variables have been considered suitable for imputation according to the defined 0.34 threshold.
#Missing values have been imputed according to KNN method.

#Processing metabolome dataset
#Original data contained 8 zeroes and 0 NAs in  5  and  0  variables, respectively. After processing 5 values present in  5 variables have been considered suitable for imputation according to the defined 0.34 threshold.
#missForest iteration 1 in progress...done!
#  missForest iteration 2 in progress...done!
#  missForest iteration 3 in progress...done!
#  missForest iteration 4 in progress...done!
#  Missing values have been imputed according to RF method.

#Processing genexpression dataset
#Original data contained 0 zeroes and 0 NAs in  0  and  0  variables, respectively. After processing no variable has been considered suitable for imputation according to the defined 0.34 threshold.

#Processing physiology dataset
#Original data contained 0 zeroes and 0 NAs in  0  and  0  variables, respectively. After processing no variable has been considered suitable for imputation according to the defined 0.34 threshold.

#REMOVING EMPTY COLUMNS OF ALL DATASETS
#Single processor core will be used. It may take a while...

#ABUNDANCE BALANCING
#Single processor core will be used. It may take a while...
#AvgIntensity balancing will be used for proteome dataset
#AvgIntensity balancing will be used for metabolome dataset
#none balancing will be used for genexpression dataset
#none balancing will be used for physiology dataset

#Processing proteome dataset
#Processing metabolome dataset
#Processing genexpression dataset
#Processing physiology dataset


#SUMMARY: Dataset imputation and balancing
#-----------------------------------------
#  proteome, metabolome, genexpression, physiology datasets were considered


#Missing value imputation.........
#proteome, metabolome were respectively imputed employing KNN, RF methods. 2 neighbors were employed for KNN calculations 
#Minimum percentage of significant values per variable/treatment to allow imputation was  0.34 
#Run in parallel TRUE 

#Abundance balancing.........
#proteome, metabolome were respectively balanced employing AvgIntensity, AvgIntensity methods.

#Job finished!

my.transformandselect.list<-transformandselect(datalist = my.preprocessed.list,initialrow = 1,initialcolumn = 2,treatment1col = 1,treatment2col = 1,treatment = 1,transf = c("Log10","Log10","none","z"),varsel = TRUE,varselthld = 0.4,varcoef = FALSE,varcoefthld = NULL)

#TRANSFORMATION OF DATASETS
#Single processor core will be used. It may take a while...
#Log10 transformation method will be used for proteome dataset
#Log10 transformation method will be used for metabolome dataset
#none transformation method will be used for genexpression dataset
#z transformation method will be used for physiology dataset

#Processing proteome dataset
#Processing metabolome dataset
#Processing genexpression dataset
#Processing physiology dataset
#SELECTING VARIABLES BASED ON CONSISTENCY
#Single processor core will be used. It may take a while...

#Processing proteome dataset
#Variable Selection based on consistency. Variable was present at least on 6  cases, or all replicates of a treatment
#Initial Variables:  1228  Selected Variables:  1228  Removed Variables:  0
#Processing metabolome dataset
#Variable Selection based on consistency. Variable was present at least on 6  cases, or all replicates of a treatment
#Initial Variables:  118  Selected Variables:  118  Removed Variables:  0
#Processing genexpression dataset
#Variable Selection based on consistency. Variable was present at least on 6  cases, or all replicates of a treatment
#Initial Variables:  6  Selected Variables:  6  Removed Variables:  0
#Processing physiology dataset
#Variable Selection based on consistency. Variable was present at least on 6  cases, or all replicates of a treatment
#Initial Variables:  5  Selected Variables:  5  Removed Variables:  0
#done!

Univariate.list<-univariate(datalist = my.transformandselect.list,initialrow = 1,initialcolumn = 2,treatment1col = 1,treatment2col = NULL,treatment = 1,parametric = TRUE,posthoc = TRUE,FDR = TRUE,round = 5)
print(Univariate.list$meansd$proteome[1:6,1:6])
print(Univariate.list$pvalue$proteome[1:6,1:6])
