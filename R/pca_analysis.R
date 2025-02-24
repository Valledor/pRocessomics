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

#' @name pca_analysis
#' @title Principal Component Analysis
#' @description A function to perform a Principal Component Analysis on the given dataset or dataset levels
#' @usage pca_analysis(datalist, annotation = NULL, initialrow = 1, initialcolumn = 3, 
#' treatment1col = 1, treatment2col = 2, treatment = 1, omiclevel = NULL)
#' @param datalist List with different preprocessed omic levels. pRoDS class object.
#' @param annotation optional, pRoAnnot class object created with pRocessomics::importannotation 
#' function containing the descriptions of the variables within the datasets
#' @param initialrow First row of numeric data within datasets.
#' @param initialcolumn First column of numeric data within datasets.
#' @param treatment1col Column number in which treatment 1 is indicated.
#' @param treatment2col Column number in which treatment 2 is indicated (if applicable).
#' If there is only one treatment indicate the same column as treatment1col.
#' @param treatment numeric. Set 1 to split the dataset according to treatment1, 2 to split the
#' dataset according to treatment 2 and 3 to split the dataset according to the combination of both treatments.
#' @param omiclevel Vector containing the name/s of omics layers to be analyzed
#' @details This function performs a PCA os scaled and centered data
#' @return A "pcaanalysis" class list containing the main outputs of the PCA which can be plotted by using pRocessomics::pca_plot() function
#' @author Luis Valledor and Laura Lamelas
#' @seealso pca_plot


#' @importFrom stats prcomp
#' @export

pca_analysis <- function(datalist, annotation = NULL, initialrow=1, initialcolumn=3, treatment1col=1, treatment2col=2, treatment=1,omiclevel=NULL){
  options(warn=-1)
  if(!inherits(datalist, "POL")) stop("A POL object is required. Please run pre-processing steps")
  datasetnames<-names(datalist)
  if(is.null(omiclevel)) omiclevel <- datasetnames
  if(any(is.na(omiclevel))) omiclevel <- datasetnames
  #if(omiclevel == "all") omiclevel <- datsetnames
  ### FLAG
  IDENTIFIER <- NULL
  
  #### Treatment table ####
  treatments <- buildtreatments(datalist,c(initialrow, initialcolumn, treatment1col, treatment2col, treatment))
  
  #### Selection of omic levels and extraction of meaningful rows and cols ####
  # Check that datasetnames exist.
  #if(any(omiclevel%in% c("c","d") ==FALSE )){stop()}
  if(any(omiclevel %in% datasetnames==FALSE)) {
    stop(
      paste("Error: Please select a valid datasetname (",
            paste(datasetnames,collapse = " "),
      ")",
      sep=""))}
  #we remove not needed rows and cols    
  filtdatalist <- lapply(datalist, function(x) x[initialrow:nrow(x),initialcolumn:ncol(x)])    
  #we generate a data matrix with selected datasets
  if(length(omiclevel)>1){  
    variablesforpca <-do.call(cbind,filtdatalist[omiclevel])
    variablenames <- lapply(filtdatalist[omiclevel], function(x) colnames(x))
    variablenames <-unlist(variablenames)
    colnames(variablesforpca) <- variablenames
    variablesforpca<-RemoveEmptyColumns(variablesforpca,1,1) #This avoids strange behaviors
  }
  if(length(omiclevel)==1){
    variablesforpca <-filtdatalist[[omiclevel]]
    variablenames <-colnames(variablesforpca)
    variablesforpca<-RemoveEmptyColumns(variablesforpca,1,1) #This avoids strange behaviors
  }
  
  #### Finishing annotations ####
  if(!is.null(annotation)){
    #### Preparation of annotations and groups ####
    annotation <- importannotation_grep(annotation)
    annotation<-subset(annotation, IDENTIFIER %in% variablenames, select=colnames(annotation))
  }
  
  #### PCA Analysis itself ####
  texts_wizard("\n\nCALCULATING PCA\n\nSingle thread computing. It may take a while...\n")
  pca_analysis <- stats::prcomp(variablesforpca,center=T,scale=T)
  
  #### Output object ####
  if(is.null(annotation)){
    results <- list(pca_analysis,treatments,names(datalist),variablesforpca)
    names(results) <- c("pca","treatments","datasetnames","originaldata")
    class(results) <- "pcaanalysis"
    texts_wizard("\npca_analysis: done!")
    return(results)
  } else {
    results <- list(pca_analysis,treatments,names(datalist),annotation,variablesforpca)
    names(results) <- c("pca","treatments","datasetnames","annotations","originaldata")
    class(results) <- "pcaanalysis"
    texts_wizard("\npca_analysis: done!")
    return(results)
  }
  
}