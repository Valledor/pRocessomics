#' @name ica_analysis
#' @title A function o perform ica analysis
#' @description This function performs the Independent Component Analysis on omic datasets
#' @usage ica_analysis(datalist, annotation = NULL, initialrow = 1, 
#' initialcolumn = 3,  treatment1col = 1, treatment2col = 2, treatment,
#' ncomp = 3, omiclevel = "all")
#' @param datalist List with different preprocessed omic levels. POL class object.
#' @param initialrow First row of numeric data within datasets.
#' @param initialcolumn First column of numeric data within datasets.
#' @param treatment1col Column number in which treatment 1 is indicated.
#' @param treatment2col Column number in which treatment 2 is indicated (if applicable).
#' If there is only one treatment indicate the same column as treatment1col.
#' @param treatment numeric. Set 1 to split the dataset according to treatment1, 2 to split the
#' dataset according to treatment 2 and 3 to split the dataset according to the combination of both treatments.
#' @param annotation optional, pRoAnnot class object created with pRocessomics::importannotation function
#' @param ncomp Integer, number of components of the model
#' @param omiclevel Vector containing the name/s of omics layers to be analyzed
#' @return icaanalysis class object 
#' @author Luis Valledor and Laura Lamelas
#' @export

#' @importFrom dplyr left_join
#' @importFrom ica icafast
 

ica_analysis <- function(datalist, annotation=NULL, initialrow=1, initialcolumn=3,  treatment1col=1, treatment2col=2,  treatment, ncomp=3, omiclevel="all"){
  options(warn=-1)
  if(!inherits(datalist, "POL")) stop("A POL object is required. Please run pre-processing steps")
  datasetnames<-names(datalist)
  if(is.null(omiclevel)) omiclevel<-datasetnames
  if(is.na(omiclevel)) omiclevel<-datasetnames
  #FLAG
  IDENTIFIER <- NULL
  #### Treatment table ####
  treatments <- buildtreatments(datalist,c(initialrow, initialcolumn, treatment1col, treatment2col, treatment))
  
  #### Selection of omic levels and extraction of meaningful rows and cols ####
  # Check that datasetnames exist.
  sapply(omiclevel, function(x) if(x %in% datasetnames==FALSE) stop(paste("Error: Please select a valid datasetname (",paste(datasetnames,collapse = " ")),")",sep=""))
  #we remove not needed rows and cols    
  filtdatalist <- lapply(datalist, function(x) x[initialrow:nrow(x),initialcolumn:ncol(x)])    
  #we generate a data matrix with selected datasets
  if(length(omiclevel)>1){  
    variablesforica <-do.call(cbind,filtdatalist[omiclevel])
    variablenames <- lapply(filtdatalist[omiclevel], function(x) colnames(x))
    variablenames <-unlist(variablenames)
    colnames(variablesforica) <- variablenames
    variablesforica<-RemoveEmptyColumns(variablesforica,1,1) #This avoids strange behaviors
  }
  if(length(omiclevel)==1){
    variablesforica <-filtdatalist[[omiclevel]]
    variablenames <-colnames(variablesforica)
    variablesforica<-RemoveEmptyColumns(variablesforica,1,1) #This avoids strange behaviors
  }
  
  
  #### Finishing annotations ####
  if(!is.null(annotation)){
    #### Preparation of annotations and groups ####
    annotation <- importannotation_grep(annotation)
    annotation<-subset(annotation, IDENTIFIER %in% variablenames, select=colnames(annotation))
  }
  
  #### ICA Analysis  itself ####
  texts_wizard("\n\nCALCULATING ICA\n\nSingle thread computing. It may take a while...\n")
  ica_analysis <- ica::icafast(variablesforica,nc=ncomp,alg="def",center=T,fun="logcosh")
  colnames(ica_analysis$W)<-variablenames
  rownames(ica_analysis$W)<-paste0("IC",1:nrow(ica_analysis$W))
  colnames(ica_analysis$S)<-paste0("IC",1:ncol(ica_analysis$S))  
  
  
  #### Output object ####
  if(is.null(annotation)){
    results <- list(ica_analysis,treatments,names(datalist),variablesforica)
    names(results) <- c("ica","treatments","datasetnames","originaldata")
    class(results) <- "icaanalysis"
    texts_wizard("\nica_analysis: done!")
    return(results)
  } else {
    results <- list(ica_analysis,treatments,names(datalist),annotation,variablesforica)
    names(results) <- c("ica","treatments","datasetnames","annotations","originaldata")
    class(results) <- "icaanalysis"
    texts_wizard("\nica_analysis: done!")
    return(results)
  }
  
}