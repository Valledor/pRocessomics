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

#' @name bda_analysis
#' @title bda analysis
#' @description A function to perform bda (block discriminant analysis)
#' @usage bda_analysis(datalist, omiclevel = NULL, annotation = NULL, 
#' initialrow = 1, initialcolumn = 3, treatment1col = 1, treatment2col = 2, 
#' treatment = 1, ncomponents = 2, keepX = NULL, folds = 3, nrepeat = 5,
#' autotune = FALSE, strength = 0.1)
#' @param datalist A "POL" object, see preprocess_omic_list for further details
#' @param omiclevel Character vector containing the name/s of omics layers to be analyzed, if set to "all" note all 
#' datasets will be concatenated and will conform the input data
#' @param annotation optional, annotations file class pRoAnnot 
#' @param initialrow First row of numeric data within datasets.
#' @param initialcolumn First column of numeric data within datasets.
#' @param treatment1col Column number in which treatment 1 is indicated.
#' @param treatment2col Column number in which treatment 2 is indicated (if applicable).
#' If there is only one treatment indicate the same column as treatment1col.
#' @param treatment numeric. Set 1 to split the dataset according to treatment1, 2 to split the
#' dataset according to treatment 2 and 3 to split the dataset according to the combination of both treatments.
#' @param ncomponents Number of components to build uop the model
#' @param keepX Number of variables to take into account for the analysis, if set to NULL 
#' all variables will be considered
#' @param folds Numeric, the folds in the Mfold cross-validation, maximum allowed 50
#' @param nrepeat Number of repeats to test the model
#' @param autotune logical, if set to TRUE several KeepX values will be tested and the better (less model error) 
#' will be emplyed
#' @param strength Numeric, accuracy of design matrix
#' @return bdaanalysis class object containing the bda analysis of input dataset
#' @author Luis Valledor and Laura Lamelas
#' @export
#' @importFrom parallel detectCores
#' @importFrom mixOmics tune.block.splsda block.splsda

bda_analysis <- function(datalist, omiclevel=NULL,annotation=NULL, initialrow=1, initialcolumn=3, treatment1col=1, treatment2col=2, 
                             treatment=1, ncomponents=2, keepX=NULL,folds=3,nrepeat=5,autotune=FALSE,strength=0.1){
  
  options(warn=-1)
  #FLAGS ----
  IDENTIFIER <- NULL
 
  if(!inherits(datalist, "POL")) stop("A POL object is required. Please run pre-processing steps")
  datasetnames<-names(datalist)
  if(is.null(omiclevel)) omiclevel <- datasetnames
  if(is.na(omiclevel)) omiclevel <- datasetnames
  if(omiclevel == "all") omiclevel <- datasetnames
  if(length(omiclevel)==1) stop("Fatal Error: You need to select at least two omic levels")
  
  #### Treatment table ####
  treatments <- buildtreatments(datalist,c(initialrow, initialcolumn, treatment1col, treatment2col, treatment))
  if(treatment==1) vtreatments <- as.vector(treatments[,1])
  if(treatment==2) vtreatments <- as.vector(treatments[,2])
  if(treatment==3) vtreatments <- as.vector(treatments[,3])
  
  #### Selection of omic levels and extraction of meaningful rows and cols ####
  # Check that datasetnames exist.
  sapply(omiclevel, function(x) if(x %in% datasetnames==FALSE) stop(paste("Error: Please select a valid datasetname (",paste(datasetnames,collapse = " ")),")",sep=""))
  datalist <- datalist[omiclevel]
  #We generate a list for analysis, remove not needed rows and cols    
  filtdatalist <- lapply(datalist, function(x) x[initialrow:nrow(x),initialcolumn:ncol(x)])    
  filtdatalist <- lapply(filtdatalist, function(x) RemoveEmptyColumns(x,1,1)) 
  datalistfordab <- filtdatalist[omiclevel] #now is redundant but I don`t change a line that works 
  variablenames <- lapply(datalistfordab, function(x) colnames(x))
  variablenames <- unlist(variablenames)
  dimlist <- unlist(lapply(lapply(filtdatalist,dim),function(x) x[1]))
  if(length(unique(dimlist))>1) stop("Cases in all datasets should be the same. Please re-check nrows")
  
  #Matrix design
  matrixdesign <- matrix(strength, ncol=length(filtdatalist), nrow =length(filtdatalist), dimnames = list(omiclevel,omiclevel))
  diag(matrixdesign) <-0
  
  #### Finishing annotations ####
  if(!is.null(annotation)){
    #### Preparation of annotations and groups ####
    annotation <- importannotation_grep(annotation)
    annotation <- subset(annotation, IDENTIFIER %in% variablenames, select=colnames(annotation))
    
  }
  
  ####Autotune####
  if(autotune==TRUE){
    #Two components is ok for most uses/users.
    testX <-lapply(datalistfordab,series)
    cat("\nTuning Block.SPLS.DA parameters. It may take a while\n")
    
    cores<-parallel::detectCores()
    set.seed(12820)
    tune.diablo <- mixOmics::tune.block.splsda(X=datalistfordab,
                                               Y=vtreatments,
                                               ncomp=ncomponents,
                                               test.keepX = testX, 
                                               design = matrixdesign, 
                                               validation = "Mfold", 
                                               folds=folds, 
                                               nrepeat=nrepeat, 
                                               dist="centroids.dist",
                                               cpus=cores)
    
    #Creamos la lista de variables a mantener
    keepX <- tune.diablo$choice.keepX
    dog <- paste("\nBlock SPLS DA autotune:\n", "The lower errors were estimated when selecting",paste(keepX,collapse=" and "), "elements (respectively for each component).",sep=" ")
    texts_wizard(dog)
    readline(prompt="Press [enter] to continue.")   
  }
  
  #### BDA Analysis and performance estimation #### 
  texts_wizard("\n\nCALCULATING BLOCK-DA\n\nSingle thread computing. It may take a while...\n")
  set.seed(12820)
  capture.output(bda_analysis <- mixOmics::block.splsda(X=datalistfordab,Y=vtreatments,ncomp=ncomponents,design = matrixdesign,keepX = keepX))

   
  #### Output object ####  
  if(is.null(annotation)){
    results <- list(bda_analysis,treatments,names(datalist),bda_analysis,datalistfordab)
    names(results) <- c("bda","treatments","datasetnames","bdanetwork","originaldata")
    class(results) <- "bdaanalysis"
    texts_wizard("\nbda_analysis: done!")
    return(results)
  } else {
    #### incluir las anotaciones  
    varsnames <- do.call(cbind,datalistfordab)
    colnames(varsnames) <- make.names(colnames(varsnames),unique = TRUE)
    set.seed(5881)
    bda_network<-mixOmics::block.splsda(X=datalistfordab,Y=vtreatments,ncomp=ncomponents,design = matrixdesign,keepX = keepX)
    results <- list(bda_analysis,treatments,names(datalist),annotation, bda_network,filtdatalist)
    names(results) <- c("bda","treatments","datasetnames","annotations","bdanetwork","originaldata")
    class(results) <- "bdaanalysis"
    texts_wizard("\nbda_analysis: done!")
    
    return(results)
  }
  
} 