#
# Copyright 2019 Luis Valledor / Laura Lamelas
# Last revision 01.2023
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

#' @name preprocess_omic_list
#' @title Preprocess datalist
#' @description A function to preprocess omics data by imputing (or not) missing values, balance abundance among samples and pre-filter the data
#' @usage preprocess_omic_list(datalist, initialrow=1, initialcolumn=2, 
#' treatment1col=1, treatment2col=2, treatment=1, imputation=NULL, imputhld=0.25, 
#' k=3, parallel=FALSE, varsel=FALSE, varselthld=0.25, abdbal=NULL)
#' @param datalist List with different preprocessed omic levels. pRoDS class object.
#' @param initialrow First row of numeric data within datasets.
#' @param initialcolumn First column of numeric data within datasets.
#' @param treatment1col Column number in which treatment 1 is indicated.
#' @param treatment2col Column number in which treatment 2 is indicated (if applicable).
#' If there is only one treatment indicate the same column as treatment1col.
#' @param treatment numeric. Set 1 to split the dataset according to treatment1, 2 to split the
#' dataset according to treatment 2 and 3 to split the dataset according to the combination of both treatments.
#' @param imputation Method for imputation, "RF" - Random Forests,"KNN" - K-nearest neighbor, or "none". If NULL Random Forest methods will be employed.This can be provided as a single value or a vector with the same number of elements than list. Each vector element can be a different imputation method or none
#' @param imputhld Maximum number of NAs, in %, per variable to be imputed.
#' @param k Number of neighbors employed in KNN algorithm.
#' @param parallel boleean. If parallel=TRUE missing value imputation algorithm will be run in parallel.
#' @param varsel boleean, If TRUE Variable selection based on consistency criteria will be used to pre-filter the data. Variables present in less than the stablished varselthld will be dropped out the analysis.
#' @param varselthld Minimum number of significant (not zero) values, in %, per variable to be kept.
#' @param abdbal Abundance balancing normalization. "sample" - sample centric approach, "AvgIntensity" - data will be processed according to a sample centric approach and then each value is multipled by average intensity (sum of the intensities of all variables within a sample) of all samples. "TreatAvgIntensity" - data is sample centric normalized and then multiplied by the average intensity of samples of specific treatments. "none" - no abundance balancing is performed. If NULL "AvgIntensity" will be employed. This can be provided as a single value or a vector with the same number of elements than list. Each vector element can be a different abundance balancing method or none
#' @return A POL class object, a list with the processed dataset 

#' @details The objective \code{preprocess_omic_list} is providing a first step of data preprocessing before uni and multivariate statistics are performed. The need of removing inconsistent variables (those present in only one or two samples, or really close to detection limit) and also defining which values are NA or zeroes for later imputation is a constant when working with omic datasets. To this end within this function the user can select the range of data to be analyzed, how it will be processed (considering all dataset together, or splitting it in the different treatments) towards the definition and imputation of missing values and also for balancing the abundance of each sample.
#' All tables must have the same structure with cases in rows and variables in columns. Row names, and columns in which treatments are defined should have the same name and order across the different datasets. By default, first columns should be devoted to indicate treatments. Row names should be unique, and preferable equal to database accession names in order to use later annotation steps (see example datasets).
#' First step of this function is deciding which value should be a 0 or NA. To this end all NAs are turned to 0, and then it is decided if specific values should remain as 0 or set as NA for imputation. \code{threshold} argument defines the maximum number of zeroes allowed for imputation within each variable and treatment (defined in \code{treatment}).
#' The second step of the function balance the abundance of the variables according to different methods: sample centric approach; average intensity of samples, or average intensity within treatment.
#' @author Luis Valledor and Laura Lamelas


#' @export

#' @importFrom methods hasArg
#' @importFrom parallel detectCores
#' @importFrom stats var



preprocess_omic_list <- function(datalist, initialrow=1, initialcolumn=2, treatment1col=1, treatment2col=2, treatment=1, imputation=NULL, imputhld=0.25, k=3, parallel=FALSE, varsel=FALSE, varselthld=0.25, abdbal=NULL){
  #### Initial checks ####
  # Dataset names, number of rows, cases, etc.
  if(methods::hasArg(datalist)==FALSE) stop("\nPlease introduce a valid dataset\n")
  #if(class(datalist) != "list")
  if(is.list(datalist)==F)
    stop("A list of matrices or dataframes corresponding to each level is expected")
  if(is.null(names(datalist)))
    stop("Names of elements in list are missing. Please provide names for each level")
  if(all(sapply(lapply(datalist, function(x) t(x[initialrow:nrow(x),initialcolumn:ncol(x)])), dim)[2,]==sapply(lapply(datalist, function(x) t(x[initialrow:nrow(x),initialcolumn:ncol(x)])), dim)[2,1])==FALSE)
    stop("The different matrices have an unequal number of individuals. The number of cases of each level should be the same. Please read instructions for more information")
  dataframesnames<-lapply(datalist, function(x) rownames(x))
  if(all(unlist(lapply(dataframesnames, function(x) identical(dataframesnames[[1]],x)))==TRUE)==FALSE)
    stop("Individuals have different names across levels. Please check original matrices")
  
  # Are all variables to be analyzed numeric?
  for(i in 1:length(datalist)){
    if(all(apply(datalist[[i]][initialrow:nrow(datalist[[i]]),initialcolumn:ncol(datalist[[i]])],2,function(x) is.numeric(x))==F)) stop(paste("Check your input. Non numeric values in",names(datalist)[i], "dataset",sep=" "))
  }
  
  # Set imputation and abundance balancing defaults if not given by user
  if(is.null(imputation)) imputation <-"RF"
  if(is.null(abdbal)) abdbal <-"AvgIntensity"
  
  # Check imputation method(s)
  if(FALSE %in% (imputation %in% c("RF","KNN","none"))==TRUE) stop("Please select a valid imputation method")
  if(length(imputation)==1) imputation <- c(rep(imputation, length(datalist)))
  if(length(imputation)!=length(datalist)) stop("A vector containing one-common- or n elements (where n is the number of datasets) indicating imputation method(s) is expected")
  
  # Check abundance balancing methods
  if(FALSE %in% (abdbal %in% c("Sample","AvgIntensity","TreatAvgIntensity","none"))==TRUE) stop("Please select a valid abundance balancing method")
  if(length(abdbal)==1) abdbal <- rep(abdbal, length(datalist))
  if(length(abdbal)!=length(datalist)) stop("A vector containing one-common- or n elements (where n is the number of datasets) indicating abundance balancing method(s) is expected")
  
  # Check k and threshold
  if(!is.numeric(k)) stop("Number of neighbors, k,  should be a numeric constant")
  if(any(imputation %in% "KNN")==FALSE) cat("Random Forest or none imputation method has been selected. k value,if provided, will be ignored") 
  threshold<-as.numeric(imputhld)
  if(!is.numeric(threshold)) stop("Threshold for defining NA or 0, threshold,  should be a numeric constant betwen 0 and 1")
  
  # Check variable selection settings
  if(varsel==TRUE & is.null(varselthld)) stop("Please select adequate threshold for VarSelect")
  if(!is.null(varselthld)) if(length(varselthld)!=1&length(varselthld)!=length(datalist)) stop("A vector containing one -common- or n (where n is the number of datasets) thresholds  is expected")
  if(!is.null(varselthld)) if(length(varselthld)==1) varselthld <- rep(varselthld, length(datalist))
  
  
  
  datasetnames<-names(datalist)

    #### Removing empty columns ####
  texts_wizard(paste("\nREMOVING EMPTY COLUMNS OF ALL DATASETS\n"))
  texts_wizard("Single processor core will be used. It may take a while...\n")
  datalist<-lapply(datalist, function(x) RemoveEmptyColumns(x,initialrow=initialrow,initialcolumn=initialcolumn))
  names(datalist)<-datasetnames
  
  
  #### Imputation ####
  datasetnames<-names(datalist)
  
  if(parallel==FALSE){
    if(!is.null(imputation)){
      texts_wizard("\n\nMISSING VALUE IMPUTATION\n")
      texts_wizard("Single processor core will be used. It may take a while...\n")
      for (i in 1:length(datalist)){
        texts_wizard(paste(imputation[i], " imputation method will be used for ", names(datalist)[i]," dataset\n",sep=""))
      }
      datalist<-mapply(function(x,y,z){
        NAorZero_Imput(x,initialrow=initialrow,initialcolumn=initialcolumn,treatment1col=treatment1col, treatment2col=treatment2col, treatment=treatment,threshold=threshold,imputation=y,k=k,cores=1,datasetname = z)
      },datalist,imputation,datasetnames,SIMPLIFY = FALSE)
    }
  }
  
  if(parallel==TRUE){
    #suppressMessages(suppressWarnings(require(doParallel)))
    cores=(parallel::detectCores()-1)
    if(!is.null(imputation)){
      texts_wizard("\n\nMISSING VALUE IMPUTATION\n")
      texts_wizard(paste("Multiple processor cores (", cores, ") will be used. It may take a while...\n",sep=""))
      for (i in 1:length(datalist)){
        texts_wizard(paste(imputation[i], " imputation method will be used for ", names(datalist)[i]," dataset\n",sep=""))
      }
        datalist<-mapply(function(x,y,z){
        NAorZero_Imput(x,initialrow=initialrow,initialcolumn=initialcolumn,treatment1col=treatment1col, treatment2col=treatment2col, treatment=treatment,threshold=threshold,imputation=y,k=k,cores=cores,datasetname = z)
      },datalist,imputation,datasetnames,SIMPLIFY = FALSE)
    }
  }
  #### Variable selection ####
  if(varsel==TRUE){
    texts_wizard("\n\nSELECTING VARIABLES BASED ON CONSISTENCY\n")
    texts_wizard("Single processor core will be used. It may take a while...\n")
    datalist<-mapply(function(x,y,z){
      VarSelect(x,initialrow=initialrow,initialcolumn=initialcolumn,treatment1col = treatment1col, treatment2col = treatment2col, treatment = treatment,datasetname = z,threshold=y)
    },datalist,varselthld,datasetnames,SIMPLIFY = FALSE)
  }
  
  #### Abundance balancing ####
  if(!is.null(abdbal)){
    texts_wizard("\n\nABUNDANCE BALANCING\n")
    texts_wizard("Single processor core will be used. It may take a while...\n")
    for (i in 1:length(datalist)){
      texts_wizard(paste(abdbal[i], " balancing will be used for ",  names(datalist)[i]," dataset\n",sep=""))
    }
    datalist<-mapply(function(x,y,z){
      AbdBal(x,initialrow=initialrow,initialcolumn=initialcolumn,treatment1col=treatment1col, treatment2col=treatment2col, treatment=treatment,norm=y,datasetname=z)
    },datalist,abdbal,datasetnames,SIMPLIFY = FALSE)
  }
  #### Final report ####
  header <- c("\n\n\nSUMMARY: Dataset imputation, filtering and balancing\n-----------------------------------------\n")
  header <- paste(header, paste(names(datalist),collapse=", "), " datasets were considered\n",sep="")
  
  class(datalist) <-"POL"
  texts_wizard("\n\n Job finished!\n\n")
  return(datalist)
}
