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
#' @name spls_analysis
#' @title sparse Partial Least Squares Analysis
#' @description This function performs a sPLS analysis on the given two levels of the dataset
#' @usage spls_analysis(datalist, annotation=NULL, initialrow=1, initialcolumn=3, 
#' treatment1col=1, treatment2col=2, treatment=1, omiclevelX, omiclevelY, keepX=NULL, 
#' keepY=NULL, autotune=FALSE, performance=FALSE, folds=3)
#' @param datalist List with different preprocessed omic levels. POL class object.
#' @param annotation optional, pRoAnnot class object created with pRocessomics::importannotation 
#' function containing the descriptions of the variables within the datasets
#' @param initialrow First row of numeric data within datasets.
#' @param initialcolumn First column of numeric data within datasets.
#' @param treatment1col Column number in which treatment 1 is indicated.
#' @param treatment2col Column number in which treatment 2 is indicated (if applicable).
#' If there is only one treatment indicate the same column as treatment1col.
#' @param treatment numeric. Set 1 to split the dataset according to treatment1, 2 to split the
#' dataset according to treatment 2 and 3 to split the dataset according to the combination of both treatments.
#' @param omiclevelX Character, quoted name of the predictor omic level (X)
#' @param omiclevelY Character, quoted name of the predicted omic level (Y)
#' @param keepX Number of variables of predictor level (X) to take into account for the analysis, if set to NULL 
#' all variables will be considered
#' @param keepY Number of variables of predicted level (Y) to take into account for the analysis, if set to NULL 
#' all variables will be considered
#' @param autotune logical, if set to TRUE several KeepX values will be tested and the better (less model error) 
#' will be employed
#' @param performance logical, analysis of the generated model performance in order to determine the 
#' optimum number of caomponents
#' @param folds numeric, the folds in the Mfold cross-validation, maximum allowed 50
#' @details The aim of this function is to perform an unsupervised analysis by predicting the one omic level (Y) using
#' quantification profiles of the predictor omic level (X). 
#' @return spls_analysis returns an object of class splsanalysis, a list that contains the following components:
#' splso: spls analysis output, 
#' treatments: considered treatments in the analysis
#' datasetnames: the names of the input dataset
#' annotations: annotations (if provided) of the elements in the analysis 
#' originaldata: original data 

#' @author Luis Valledor and Laura Lamelas
#' @seealso da_analysis, da_plot, spls_plot
#' 


#' @importFrom mixOmics spls perf
#' @importFrom utils txtProgressBar
#' @importFrom plyr count
#' @export

spls_analysis <- function(datalist, annotation=NULL, initialrow=1, initialcolumn=3, treatment1col=1, treatment2col=2, 
                          treatment=1, omiclevelX, omiclevelY, keepX=NULL, keepY=NULL,autotune=FALSE, performance=FALSE, folds=3){
  options(warn=-1)
  
  IDENTIFIER <- NULL
  
  if(!inherits(datalist, "POL")) stop("A POL object is required. Please run pre-processing steps")
  datasetnames<-names(datalist)
  if(is.null(omiclevelX)) stop("Select omic level for prediction matrix (omiclevelX)")
  if(is.null(omiclevelY)) stop("Select omic level for response matrix (omiclevelY)")
  sapply(omiclevelX, function(x) if(x %in% datasetnames==FALSE) stop("Error: Please select a valid datasetname for prediction matrix"))
  sapply(omiclevelY, function(x) if(x %in% datasetnames==FALSE) stop("Error: Please select a valid datasetname for response matrix"))
  
  
  #### Treatment table ####
  treatments <- buildtreatments(datalist,c(initialrow, initialcolumn, treatment1col, treatment2col, treatment))
  if(treatment==1) vtreatments <- as.vector(treatments[,1])
  if(treatment==2) vtreatments <- as.vector(treatments[,2])
  if(treatment==3) vtreatments <- as.vector(treatments[,3])
  
  #### Building of predictor and response matrices ####
  
  filtdatalist <- lapply(datalist, function(x) x[initialrow:nrow(x),initialcolumn:ncol(x)])    
  #we generate a data matrix with selected datasets
  if(length(omiclevelX)>1){  
    predictormatrix <-do.call(cbind,filtdatalist[omiclevelX])
    variablenames <- lapply(filtdatalist[omiclevelX], function(x) colnames(x))
    variablenames <-unlist(variablenames)
    colnames(predictormatrix) <- variablenames
    predictormatrix<-RemoveEmptyColumns(predictormatrix,1,1) #This avoids strange behaviors
  }
  if(length(omiclevelX)==1){
    predictormatrix <-filtdatalist[[omiclevelX]]
    variablenames <-colnames(predictormatrix)
    predictormatrix<-RemoveEmptyColumns(predictormatrix,1,1) #This avoids strange behaviors
  }
  
  if(length(omiclevelY)>1){  
    responsematrix <-do.call(cbind,filtdatalist[omiclevelY])
    variablenamesY <- lapply(filtdatalist[omiclevelY], function(x) colnames(x))
    variablenamesY <-unlist(variablenamesY)
    colnames(responsematrix) <- variablenamesY
    responsematrix<-RemoveEmptyColumns(responsematrix,1,1) #This avoids strange behaviors
  }
  if(length(omiclevelY)==1){
    responsematrix <-filtdatalist[[omiclevelY]]
    variablenamesY <-colnames(responsematrix)
    responsematrix<-RemoveEmptyColumns(responsematrix,1,1) #This avoids strange behaviors
  }
  
  originaldata <- cbind(predictormatrix,responsematrix)
  
  #### Building annotation table #### 
  if(!is.null(annotation)){
    #### Preparation of annotations and groups ####
    annotation <- importannotation_grep(annotation)
    annotationX<-subset(annotation, IDENTIFIER %in% variablenames, select=colnames(annotation))
    annotationY<-subset(annotation, IDENTIFIER %in% variablenamesY, select=colnames(annotation))
    variablenamesX <- make.names(strtrim(annotationX[,2],50),unique=TRUE)
    variablenamesY <- make.names(strtrim(annotationY[,2],50),unique=TRUE)
 
  }
  
  #### Autotune SPLS ####
  if(autotune==TRUE){
    miperf <- c()
    vvtX <-series(predictormatrix)
    vvtY <-series(responsematrix)
    texts_wizard("\n\nAutotune in progress.\n")
    
    pb <- utils::txtProgressBar(min = 0, max = length(vvtX), style = 3)
    dog <-paste("\nTesting different sPLS models. Keeping ", paste(vvtX,collapse=", ")," and ", paste(vvtY,collapse=", "), " most important variables for predictor and response matrices will be tested.\n Single core will be used so it may take a while...")
    
    for(i in 1:length(vvtX)){
      for(j in 1:length(vvtY)){
        set.seed(5881)
        spls_tune <- mixOmics::spls(predictormatrix,responsematrix,mode ="regression",ncomp=2,keepX=c(vvtX[i],vvtX[i]),keepY = c(vvtY[j],vvtY[j]))
        set.seed(5881)
        utils::capture.output(perf <- mixOmics::perf(spls_tune, validation="Mfold", folds=3, criterion="Mfold", progressbar=FALSE,nrepeat = 50))
        miperf <- c(miperf,sum(perf$Q2.total))
      }
      utils::setTxtProgressBar(pb, i)
    }
    
    perfmatrix <- matrix(miperf,ncol = length(vvtY),byrow = T)
    row.names(perfmatrix)<-vvtX
    colnames(perfmatrix)<-vvtY
    
    selection<- as.vector(which(perfmatrix == min(perfmatrix), arr.ind = TRUE))
    selection<-c(vvtX[selection[1]],vvtY[selection[2]])
    keepX <- c(selection[1],selection[1])
    keepY <- c(selection[2],selection[2])
    
    dog <- paste("\nLowest error (total Q2):",min(miperf), "Was obtained when selecting: ", selection[1], " and ", selection[2],  "variables of prediction and response matrix respectively")
    texts_wizard(dog)
    readline(prompt="Press [enter] to continue")
    
  }
  
  #### DA Analysis and performance estimation #### 
  texts_wizard("\nCALCULATING sPLS\nIt may take a while...\n")
  set.seed(12820)
  spls_analysis <- mixOmics::spls(predictormatrix,responsematrix,mode ="regression",ncomp=2,keepX=keepX,keepY =keepY)
  
  if(performance==TRUE){      
    texts_wizard("\nAnalyzing model performance\nIt may take a while...\n")
    tentfolds <- min(plyr::count(as.vector(vtreatments))[,2]) #Sacamos el fold maximo para que no pete performance
    if(is.null(folds)) folds<-tentfolds
    if(folds>tentfolds) folds<-tentfolds
    set.seed(12820)
    perform <- mixOmics::perf(spls_analysis, validation="Mfold", folds=folds, criterion="Mfold", progressbar=TRUE,nrepeat = 50)
    dog <- paste("\nQ2 Total errors: PC1:",perform$Q2.total[1],"PC2:",perform$Q2.total[2],sep=" ")
    texts_wizard(dog)
  }
  #### Output object ####
  # without annotations
  if(is.null(annotation)){
    resultado <- list(spls_analysis,treatments,names(datalist),originaldata)
    names(resultado) <- c("splso","treatments","datasetnames","originaldata")
    class(resultado) <- "splsanalysis"
    cat("\nDone!")
    return(resultado)
  } else {
    # Annotations. This is horrible, but makes the trick to get on screen networks with names
    #colnames(predictormatrix) <-variablenamesX
    #colnames(responsematrix) <- variablenamesY
    set.seed(12820)
    spls_network <- mixOmics::spls(predictormatrix,responsematrix,mode ="regression",ncomp=2,keepX=keepX,keepY=keepY)
    resultado <- list(spls_analysis,treatments,names(datalist),annotation,spls_network,originaldata)
    names(resultado) <- c("splso","treatments","datasetnames","annotations","splsnetwork","originaldata")
    class(resultado) <- "splsanalysis"
    cat("\nDone!")
    return(resultado)
  }
}
