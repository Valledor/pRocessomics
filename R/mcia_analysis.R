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

#' @name mcia_analysis
#' @title Multiple coinertia analysis
#' @description A function to perform a multiple coinertia analysis
#' @usage mcia_analysis(datalist, annotation = NULL, initialrow = 1, initialcolumn = 3,
#' treatment1col = 1, treatment2col = 2, treatment = 1, omiclevel = NULL, keptaxes = 5)
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
#' @param omiclevel Vector containing the name/s of omics layers to be analyzed
#' @param keptaxes Integer, numeber of omiclevels to be kept for the analysis
#' @author Luis Valledor and Laura Lamelas
#' @return A "mcoaanalysis" class list object containing the result of the mcia analysis 
#' @importFrom omicade4 mcia
#' @export



mcia_analysis <- function(datalist, annotation=NULL, initialrow=1, initialcolumn=3,treatment1col=1, treatment2col=2, treatment=1, omiclevel=NULL, keptaxes=5){
  options(warn=-1)
  #FLAG
  IDENTIFIER <- NULL
  
  if(!inherits(datalist, "POL")) stop("A POL object is required. Please run pre-processing steps")
  datasetnames<-names(datalist)
  if(is.null(omiclevel)) omiclevel<-datasetnames
  #if(is.na(omiclevel)) omiclevel<-datasetnames
  if(omiclevel == "all") omiclevel <- datasetnames
  if(length(omiclevel)==1) stop("Fatal Error: You need to select at least two omic levels for performing this analysis")
  if (treatment %in% c(1,2,3) == FALSE) stop("Treatment should be 1, 2 or 3")
  #### Treatment table ####
  treatments <- buildtreatments(datalist,c(initialrow, initialcolumn, treatment1col, treatment2col, treatment))
  
  #### Selection of omic levels and extraction of meaningful rows and cols ####
  # Check that datasetnames exist.
  sapply(omiclevel, function(x) if(x %in% datasetnames==FALSE) stop(paste("Error: Please select a valid datasetname (",paste(datasetnames,collapse = " ")),")",sep=""))
  #we remove not needed rows and cols    
  filtdatalist <- lapply(datalist, function(x) x[initialrow:nrow(x),initialcolumn:ncol(x)])    
  #We generate a list for analysis, removing not needed rows and cols    
  filtdatalist <- lapply(datalist, function(x) x[initialrow:nrow(x),initialcolumn:ncol(x)])    
  datalistformcia <-filtdatalist[omiclevel]
  datalistformcia <- lapply(datalistformcia, function(x) RemoveEmptyColumns(x,1,1)) #This avoids strange behaviors 
  variablenames <- lapply(datalistformcia, function(x) colnames(x))
  variablenames <-unlist(variablenames)
  dimlist<-unlist(lapply(lapply(filtdatalist,dim),function(x) x[1]))
  if(length(unique(dimlist))>1) stop("Cases in all datasets should be the same. Please re-check nrows")

  #### Finishing annotations ####
  if(!is.null(annotation)){
    #### Preparation of annotations and groups ####
    annotation <- importannotation_grep(annotation)
    annotation<-subset(annotation, IDENTIFIER %in% variablenames, select=colnames(annotation))
  }
  
  #*#*#* FOR MCIA ANALYSIS, CASES IN COLUMNS, VARIABLES IN ROWS
  datalistformcia<-lapply(datalistformcia, function(x) t(x))

  #### MCIA Analysis itself ####
  texts_wizard("\n\nCALCULATING MCIA\n\nSingle thread computing. It may take a while...\n")
  mcia_analysis <- omicade4::mcia(datalistformcia, cia.nf=keptaxes, cia.scan=F) 
  colnames(mcia_analysis$mcoa$Tco)<-colnames(mcia_analysis$mcoa$SynVar)
  colnames(mcia_analysis$mcoa$Tl1)<-colnames(mcia_analysis$mcoa$SynVar)
  colnames(mcia_analysis$mcoa$cov2)<-paste0("Cov2_",colnames(mcia_analysis$mcoa$SynVar))
  datalistformcia<-lapply(datalistformcia, function(x) t(x)) #we set original dataset to previous form

  #### Output object ####
  if(is.null(annotation)){
    results <- list(mcia_analysis,treatments,names(datalist),datalistformcia)
    names(results) <- c("mcia","treatments","datasetnames","originaldata")
    class(results) <- "mciaanalysis"
    texts_wizard("\nmcia_analysis: done!")
    return(results)
  } else {
    results <- list(mcia_analysis,treatments,names(datalist),annotation,datalistformcia)
    names(results) <- c("mcia","treatments","datasetnames","annotations","originaldata")
    class(results) <- "mciaanalysis"
    texts_wizard("\nmcia_analysis: done!")
    return(results)
  }
  
}