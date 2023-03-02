
#' @name da_analysis
#' @title A function to perform a sparce partial least squared discriminant analysis
#' @description This function takes a preprocessed omic datalist and performs a spls-da, multivariate discriminant analysis
#' @usage da_analysis (datalist, annotation = NULL, initialrow=1, 
#' initialcolumn=3, treatment1col=1, treatment2col=2, treatment=1, 
#' omiclevel="all",keepX=NULL, autotune=FALSE, performance=FALSE, folds=3)
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
#' @param omiclevel Vector containing the name/s of omics layers to be analyzed.
#' @param keepX Number of variables to take into account for the analysis, if set to NULL 
#' all variables will be considered
#' @param autotune logical, if set to TRUE several KeepX values will be tested and the better (less model error) 
#' will be employed
#' @param performance logical, analysis of the generated model performance in order to determine the 
#' optimum number of caomponents
#' @param folds numeric, the folds in the Mfold cross-validation, maximum allowed 50
#' @details The aim of this function is to perform a supervised analysis by combining the assayed 
#' experimental treatments against element values such as proteins/metabolites/transcripts 
#' quantification profiles. 
#' @return da_analysis returns an object of class daanalysis, a list that contains the following components:
#' da: spls-da analysis output, 
#' treatments: considered treatments in the analysis
#' datasetnames: the names of the input dataset
#' annotations: annotations (if provided) of the elements in the analysis 
#' originaldata: original data 

#' @author Luis Valledor and Laura Lamelas
#' @seealso da_plot da_analysis_wizard
#' 
#' @export
#' @importFrom utils txtProgressBar capture.output
#' @importFrom mixOmics splsda perf color.mixo
#' @importFrom plyr count
#' @importFrom graphics plot

da_analysis <- function(datalist, annotation=NULL, initialrow=1, initialcolumn=3, treatment1col=1, treatment2col=2, 
                        treatment=1, omiclevel="all", keepX=NULL, autotune=FALSE, performance=FALSE, folds=3){
  options(warn=-1)
  #FLAGS ----
  IDENTIFIER <- NULL
  
  if(!inherits(datalist, "POL")) stop("A POL object is required. Please run pre-processing steps")
  datasetnames<-names(datalist)
  if(omiclevel == "all") omiclevel <- datasetnames
  if(is.null(omiclevel)) omiclevel <- datasetnames
  if(is.na(omiclevel)) omiclevel <- datasetnames
  
  #### Treatment table ####
  treatments <- buildtreatments(datalist,c(initialrow, initialcolumn, treatment1col, treatment2col, treatment))
  if(treatment==1) vtreatments <- as.vector(treatments[,1])
  if(treatment==2) vtreatments <- as.vector(treatments[,2])
  if(treatment==3) vtreatments <- as.vector(treatments[,3])
  
  #### Selection of omic levels and extraction of meaningful rows and cols ####
  # Check that datasetnames exist.
  sapply(omiclevel, function(x) if(x %in% datasetnames==FALSE) stop(paste("Error: Please select a valid datasetname (",paste(datasetnames,collapse = " ")),")",sep=""))
  #we remove not needed rows and cols    
  filtdatalist <- lapply(datalist, function(x) x[initialrow:nrow(x),initialcolumn:ncol(x)])    
  #we generate a data matrix with selected datasets
  if(length(omiclevel)>1){  
    variablesforda <-do.call(cbind,filtdatalist[omiclevel])
    variablenames <- lapply(filtdatalist[omiclevel], function(x) colnames(x))
    variablenames <-unlist(variablenames)
    colnames(variablesforda) <- variablenames
    variablesforda<-RemoveEmptyColumns(variablesforda,1,1) #This avoids strange behaviors
  }
  if(length(omiclevel)==1){
    variablesforda <-filtdatalist[[omiclevel]]
    variablenames <-colnames(variablesforda)
    variablesforda<-RemoveEmptyColumns(variablesforda,1,1) #This avoids strange behaviors
  }
  
  
  #### Finishing annotations ####
  if(!is.null(annotation)){
    #### Preparation of annotations and groups ####
    annotation <- importannotation_grep(annotation)
    annotation<-subset(annotation, IDENTIFIER %in% variablenames, select=colnames(annotation))
      }
  
  
  ####Autotune####
  
  if(autotune==TRUE){
    miperf <- c()
    vvt <-series(variablesforda)
    texts_wizard("\n\nAutotune in progress.\n")
    pb <- utils::txtProgressBar(min = 0, max = length(vvt), style = 3)
    dog <-paste("\nTesting different DA models. Keeping ", paste(vvt,collapse=", ")," most important variables will be tested.\n Single core will be used so it may take a while...")
    for(i in 1:length(vvt)){
      set.seed(5881)
      da_tune <- mixOmics::splsda(variablesforda,as.factor(vtreatments),ncomp=2,keepX=c(vvt[i],vvt[i]))
      set.seed(5881)
      utils::capture.output(perform <- mixOmics::perf(da_tune, validation="Mfold", folds=3, criterion="Mfold", progressbar=FALSE,nrepeat = 50))
      miperf <-cbind(miperf, perform$error.rate$BER[,1])
      utils::setTxtProgressBar(pb, i)
    }
    xc1<-min(which(miperf[1,]==min(miperf[1,]))) 
    xc2<-min(which(miperf[2,]==min(miperf[2,])))
    keepX <-c(xc1,xc2)
    dog<-paste("\nDA autotune:\n", "The lower errors,",paste(round(miperf[keepX],2),collapse=" and "),", were estimated when selecting",paste(vvt[keepX],collapse=" and "), "elements for Comp1 and Comp2, respectively.",sep=" ")
    texts_wizard(dog)
    readline(prompt="Press [enter] to continue.")
  }
  
  #### DA Analysis and performance estimation #### 
  texts_wizard("\n\nCALCULATING DA\n\nSingle thread computing. It may take a while...\n")
  set.seed(5881)
  da_analysis <- mixOmics::splsda(variablesforda,as.factor(vtreatments),ncomp=2,keepX =keepX)
  texts_wizard("\nFinished!\n")
  
  if(performance==TRUE){
    cat("\nAnalyzing model performance\nIt may take a while...\n")
    tentfolds <- min(plyr::count(vtreatments)[,2]) #Max fold that can be used for avoiding crashing
    if(is.null(folds)) folds<-tentfolds
    if(folds>tentfolds) folds<-tentfolds
    set.seed(5881)
    perform <- mixOmics::perf(da_analysis, validation="Mfold", folds=folds, criterion="Mfold", progressbar=TRUE,nrepeat = 50)
    graphics::plot(perform, col = mixOmics::color.mixo(5:7), sd = T, legend.position = "horizontal") 
  }
  
  #### Output object ####  
  if(is.null(annotation)){
    results <- list(da_analysis,treatments,names(datalist),da_analysis,variablesforda)
    names(results) <- c("da","treatments","datasetnames","danetwork","originaldata")
    class(results) <- "daanalysis"
    texts_wizard("\nda_analysis: done!")
    return(results)
  } else {
    colnames(variablesforda) <- make.names(colnames(variablesforda),unique = TRUE)
    set.seed(5881)
    da_network <- mixOmics::splsda(variablesforda,as.factor(vtreatments),ncomp=2,keepX =keepX)
    results <- list(da_analysis,treatments,names(datalist),annotation, da_network,variablesforda)
    names(results) <- c("da","treatments","datasetnames","annotations","danetwork","originaldata")
    class(results) <- "daanalysis"
    texts_wizard("\nda_analysis: done!")
    return(results)
  }
  
}