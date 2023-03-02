#
# Copyright 2019 Luis Valledor / Laura Lamelas
# Last revision 31.03.2020
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

#' @name preprocess_wizard
#' @title Wizard for easily preprocess omic datasets
#' @description This function is used to define the parameters and preprocess omic dataset lists in a friendly way.
#' @usage preprocess_wizard(outputobjectname)
#' @param outputobjectname quoted name to assign to the preprocessed object
#' @details This function pre-processes your data by handling missing values, differences in relative sample abundance and prefiltering. All the mentioned steps are optional and the user will be able to choose during the wizard run.
#' Please check \url{https://github.com/Valledor/pRocessomics/} for more information and examples.
#' @return A POL class list containing preprocessed omics data
#' @author Luis Valledor and Laura Lamelas
#' @seealso \code{\link{preprocess_omic_list}}.     
#'        
#' @importFrom utils select.list 
#' @importFrom sjmisc is_empty
#' @export


preprocess_wizard <- function(outputobjectname=NULL){
  #### Starting console logging ####
  diahora <-gsub(":","",gsub(" ","_", gsub("-","",Sys.time())))
  diahora <-substr(diahora,3,nchar(diahora)-2)
  sinkfilename <- paste(diahora,"_PreProcWzrdLog.txt",sep="")
  sink(sinkfilename, split = T)
  
  
  ## Defining initial values
  k=3
  threshold=0.25
  
  
  #### Initial texts ####
  headings_wizards("Welcome to pRocessomics Data Preprocessing Wizard")
  texts_wizard("\n\nThis wizard is aimed to guide you for properly preprocess your datasets. Please note that this module will only work if you have imported your data using importfromexcel(). Please read 'https://github.com/Valledor/pRocessomics/wiki/Pre-processing-your-data' before continuing.\n")
  texts_wizard("If you have imported the data by yourself please change the class of your datalist to 'pRoDS' before proceeding. Since no data structure and quality check will be performed. We strongly recommend using importfromexcel().\n")
  texts_wizard("Press [enter] to continue...")
  invisible(readline())
  
  
  #### Select datalist ####
  # We ask the user to input one datalist to be used
  datalists <- Filter(function(x) "pRoDS" %in% class(get(x)), ls(envir = .GlobalEnv))
  if(length(datalists)==0) stop("No valid pRoDS objects were found. Please import your data with importfromexcel().")
  texts_wizard("\nWe have detected the following datasets available. Please enter the number of the dataset do you want to process and press [enter]. If your dataset is not listed please exit pressing [esc] and re-import your data.")
  datalists<- utils::select.list(datalists, preselect = NULL, multiple = F, title = NULL, graphics = FALSE)
  datalist <-get(datalists,envir=.GlobalEnv)
  dog<-as.character(paste("You have selected",datalists,"dataset. This dataset has the following omic levels:\n",sep=" "))
  texts_wizard(dog)
  cat(names(datalist))
  
  if(exists(paste(datalists,"_dsettings",sep=""), envir = .GlobalEnv)==TRUE){
    settings <- get(paste(datalists,"_dsettings",sep=""),envir = .GlobalEnv)
    initialrow <- settings[1]
    initialcolumn <- settings[2]
    treatment1col <- settings[3]
    treatment2col <- settings[4]
    texts_wizard("\n\nWe have detected information about previous import.")
  
    dog<-as.character(paste("\nInitialrow: ", initialrow,"\nInitialcolumn: ", initialcolumn, "\nTreatment1: ", treatment1col, "\nTreatment2:", treatment2col))
    texts_wizard(dog)
    texts_wizard("\nPress [enter] to continue...")
    invisible(readline())
    if(is.na(treatment2col)==T) treatment2col <-NULL
  } else {
    texts_wizard("\n\nInformation about dataset structure is not available. This is probably the result of a manual importation of data. However, if you employed importfromexcel wizard something went wrong and you need to repeat importing process.")
    texts_wizard("\nDo you want to continue and set this values?")
    respuestaimportacion <- utils::select.list(choices = truefalse)
    if(respuestaimportacion == "Yes") respuestaimportacion <- TRUE
    if(respuestaimportacion == "No") respuestaimportacion <- FALSE
    
    
    if(respuestaimportacion==FALSE) stop("Preprocess Wizard aborted by user")
    
    texts_wizard("\nBelow you will find a small portion of your dataset so you can easily check which is your starting row, starting column, and treatment columns.\n")
    datasetsample<-datalist[[1]][1:5,1:5]#[1:7,1:7]
    datasetsample<-rbind(colnames(datasetsample),datasetsample)
    datasetsample<-cbind(rownames(datasetsample),datasetsample)
    rownames(datasetsample) <- c(" ", "Row 1","Row 2","Row 3","Row 4","Row 5")
    colnames(datasetsample) <- c(" ", "Col 1","Col 2","Col 3","Col 4","Col 5")
    print(datasetsample)
   
    treatment1col <- as.integer(readline(prompt="Enter the number of the column with the first treatment and press [enter]: "))
    if(is.null(treatment1col)) print("null")
    if(is.na(treatment1col)) print("NA")
    if(is_empty(treatment1col)) print("empty")
    checkuserinput(treatment1col ,nrow=NULL,ncol(datalist[[1]]))
    treatment2col <-  readline(prompt="Enter the number of the column with the second treatment and press [enter]: (optional, press [enter] if you only have one treatment/level): ")
    if(treatment2col!="") {treatment2col <- as.integer(treatment2col)
    checkuserinput <-checkuserinput(treatment2col,nrow=NULL,ncol(datalist[[1]]))
    if(treatment2col==treatment1col) treatment2col=NULL}
    if(treatment2col=="") treatment2col <- NULL
    initialcolumn <- as.integer(readline(prompt="Enter the number of the column with the first variable (numerical, not treatment) and press [enter]:  "))
    checkuserinput(initialcolumn ,nrow=NULL,ncol(datalist[[1]]))
    #variablenamescolumn <- as.integer(readline(prompt="Enter the number of the row containing variable names and press [enter]: "))
    #checkuserinput(variablenamescolumn,nrow(datalist[[1]]),ncol=NULL) 
    initialrow <- as.integer(readline(prompt="Enter the number of the row with the first case (individual) and press [enter]: "))
    checkuserinput(initialrow,nrow(datalist[[1]]),ncol=NULL)
    
    #### Remove all rows but initial ####
    if(is.null(initialrow)){
      texts_wizard("\n\nWARNING: initial row was not provided. Considering initialrow=1 (Data of first case is in row 1). Please cancel and set proper initialrow if you don't like this setting.\n\n")
      readline(prompt="Press [enter] to continue")
      initialrow=1}
    if(initialrow!=1){
      texts_wizard("\n\n Removing rows without any relevant data for downstream analyses (rows 1:initialrow).\n")
      dataset <- lapply(dataset, function(x){x <- x[initialrow:nrow(x),]
      return(x)
      })
      initialrow=1
    }
  }

    if(!is.null(treatment2col)){
      texts_wizard("\nPlease select which treatment you want to use for splitting your data.\n")
    treatment<-utils::select.list(c("Treatment 1", "Treatment 2","Combination of Treatments 1 and 2"), preselect = NULL, multiple = F, title = NULL, graphics = FALSE)
    if(treatment=="Treatment 1") treatment<-1
    if(treatment=="Treatment 2") treatment<-2
    if(treatment=="Combination of Treatments 1 and 2") treatment <-3
    }
    

  #### Imputation ####
texts_wizard("\nNow we are ready to start processing the dataset. The first decision to be taken is how you will deal with missing values (NAs). NAs indicate that this variable is not present in your dataset or they are artifacts that need to be corrected? You can either change NAs to 0 (assuming that they are below detection limit, so 0 is OK) or perform imputation (a mathematical approach to assing proper figures to NAs). Please check 'https://github.com/Valledor/pRocessomics/wiki/Pre-processing-your-data' if you feel you need to know more about missing value imputation before continuing.")
texts_wizard("\nDo you want to perform missing value imputation (Yes) or transform all your missing values into 0s? (No)")
imputation <- utils::select.list(choices = truefalse)
if(imputation == "Yes") imputation <- TRUE
if(imputation == "No") imputation <- FALSE

if(imputation==F){
  imputation="none"
  treatment=1
  }


if(imputation==T) { #### Imputation==TRUE####
  texts_wizard("\nData will be split considering the treatments before imputation. Then you must choose the imputation threshold (maximum percentage of values to be imputed within each treatment). Imputation threshold below 20% (<= 0.2) is recommended, however if you have 3 to 5 replicates per treatment level you can set it to 0.34 (to impute maximum one value). Check 'https://github.com/Valledor/pRocessomics/wiki' for more info.")
  texts_wizard("\nNot imputed values will be converted to 0s.")
if(!is.null(treatment2col)) {
  texts_wizard("\nSince you have defined more than one treatment, its time to choose how will you split your data. If you don't know now what to do please 'https://github.com/Valledor/pRocessomics/wiki/Pre-processing-your-data'. It is a 5 min reading that will teach you the basics for choicing properly.")
}

if(!is.null(treatment2col)) {
  #cat("\nPlease select which treatment you want to use for splitting your data.\n")
  dog<-as.character(paste("\nData will be split according to", treatment, "This treatment has these levels:"))
  texts_wizard(dog)
  }
  
if(is.null(treatment2col)) {
  texts_wizard("\nData will be split according to treatment1 (unique treatment). This treatment has these levels: ")
  treatment<-1
}

#### Check if treatments are OK ####
# First, treatment columns should be before (left to right) data columns
if(treatment1col>=initialcolumn) stop("TABLE STRUCTURE ERROR: Treatment 1 column should be located before first data column (considering left to right order). Re-check format of your table. Please see documentation or visit wiki at https://github.com/Valledor/pRocessomics to lear how to format your datasets to be used in pRocessomics")
if(!is.null(treatment2col)) if(treatment2col>=initialcolumn) stop("TABLE STRUCTURE ERROR: Treatment 2 column should be located before first data column (considering left to right order). Re-check format of your table. Please see documentation or visit wiki at https://github.com/Valledor/pRocessomics to lear how to format your datasets to be used in pRocessomics")
# segunda asuncion, han de ser factorizables
if(treatment==1) treatmentf <- factor(datalist[[1]][initialrow:nrow(datalist[[1]]),treatment1col],levels=unique(datalist[[1]][initialrow:nrow(datalist[[1]]),treatment1col]))
if(treatment==2) treatmentf <- factor(datalist[[1]][initialrow:nrow(datalist[[1]]),treatment2col],levels=unique(datalist[[1]][initialrow:nrow(datalist[[1]]),treatment2col]))
if(treatment==3) treatmentf <- factor(paste(datalist[[1]][initialrow:nrow(datalist[[1]]),treatment1col],datalist[[1]][initialrow:nrow(datalist[[1]]),treatment2col],sep="_"), levels=unique(paste(datalist[[1]][initialrow:nrow(datalist[[1]]),treatment1col],datalist[[1]][initialrow:nrow(datalist[[1]]),treatment2col],sep="_")))

#### Recommend a threshold ####
replicatesperlevel<-length(treatmentf)/length(levels(treatmentf))
recomthreshold <- round(replicatesperlevel*0.2)/round(replicatesperlevel)
recomthreshold <- round(recomthreshold + 5*10^(-3), 2)
dog<-as.character(paste(c(levels(treatmentf)), collapse =" "))
texts_wizard(dog)
dog<-as.character(paste("\nThe average number of cases per treatment level is", round(length(treatmentf)/length(levels(treatmentf))), ", in consequence we recommend the imputation of a maximum of", round(length(treatmentf)*0.2/length(levels(treatmentf))),"value(s) which corresponds to a", recomthreshold, "threshold.\n"))
texts_wizard(dog)
threshold <- recomthreshold
thresholdd <- readline("Press [enter] to accept suggested threshold, or introduce a value [0:1] and press [enter]:")

if(thresholdd!="") {
  threshold <- as.numeric(threshold)
  if((threshold>=0&&threshold<=1)==FALSE) stop("FATAL ERROR: Threshold should be a decimal value between 0 and 1")
}

#### Select imputation method ####
texts_wizard("\nNow you need to define which imputation method do you want to use. pRocessomics provides K-Nearest-Neighbour (KNN) and Random Forests (RF). We recommend to use RF, but in doubt please check 'https://github.com/Valledor/pRocessomics/wiki' for more information.")
texts_wizard("\nIntroduce the number of the method of your choice and press [enter]:")
imputation <- utils::select.list(c("RF","KNN"), multiple = F, title = NULL, graphics = FALSE)

if(imputation=="KNN"){
  k <- readline("Introduce the number of neighbours, k, and press [enter]: ")
  k <- as.numeric(k)
  if(sjmisc::is_empty(k)==T){stop("ERROR. Please re-run this wizard and introduce and provide a valid k value.")}
  if(k%%1>0) stop("ERROR. k should be an integer. Please re-run this wizard and introduce a valid k value.")
  if(k>replicatesperlevel*2) stop("k is out of range. Introduce a smaller value of k. Hint: you should not employ a k greater than the number of replicates per treatment level you have -1.")
}
}

if(imputation==F) imputation <- "none"
k=3 # This should be here to define k

  #### VarSel ####
texts_wizard("\n\nDo you want to apply a consistency pre-filter to your data? This will discard those variables that are not present in a suffient number of samples and can introduce noise in later analyses. If you prefer to keep your dataset as it is select No.")
varsel=F
varselthld<-0.2
varsel <- utils::select.list(choices = truefalse)
if(varsel == "Yes") varsel <- TRUE
if(varsel == "No") varsel <- FALSE

if(varsel==T){
  texts_wizard("\nNow you have to establish the consistency threshold. This is the minimum proportion of significant values to satisfy consistency (e.g. threshold = 0.2 will drop out variables present in less than 20% of the cases). Note that qualitative variables (those present in all the replicates of a treatment) will NOT be discarded independently of the applied threshold.")
  varselthld<-readline("Introduce the threshold, a number between 0-1 and press [enter]: ")
  varselthld<-as.numeric(varselthld)
  if(sjmisc::is_empty(varselthld)==T){stop("ERROR. Please re-run this wizard and introduce and provide a valid threshold.")}
  if(dplyr::between(varselthld,0,1)==F){stop("ERROR. Please re-run this wizard and introduce and provide a threshold value within 0 and 1.")}
}

  #### Abdbal ####
texts_wizard("\n\nFinally it is time to decide if you want to balance your data. You can perform no-balancing (none), sample-centric (Sample), sample-average (AvgIntensity), or treatment-average intensity (TreatAvgIntensity) balancing. Usually you will use Sample (each variable value is represented as parts per unit considering total sample abundance) or AvgIntensity (Sample-centric multiplied by samples average intensity. Gives more realistic values). Please check 'https://github.com/Valledor/pRocessomics/wiki/Pre-processing-your-data' for further information. \nPlease select which balancing method you want to apply")
choices <- c("none","Sample","AvgIntensity","TreatAvgIntensity")
abdbal<- utils::select.list(choices, preselect = NULL, multiple = F, title = NULL, graphics = FALSE)

  
  #### Datasets to be analyzed ####
texts_wizard("\nSelect which omic levels you want to process with this settings:")
datalistselection <- select.list.edit(names(datalist), preselect = NULL, multiple = T, title = NULL, graphics = FALSE)
finaldatalist <- datalist[datalistselection]


  #### Parallelization ####
texts_wizard("\nDo you want to employ parallel computing to speed up processing?")
parallel <- utils::select.list(choices = truefalse)
if(parallel == "Yes") parallel <- TRUE
if(parallel == "No") parallel <- FALSE

  

  #### Launch functions ####
#return(c(initialrow, initialcolumn, treatment1col, treatment2col, treatment, imputation, abdbal, threshold, k, parallel))
myprocesseddata <- preprocess_omic_list(datalist=finaldatalist, initialrow=initialrow, initialcolumn=initialcolumn, treatment1col=treatment1col, treatment2col=treatment2col, treatment=treatment, imputation=imputation, imputhld =threshold, k=k, parallel=parallel, abdbal= abdbal, varsel=varsel, varselthld=varselthld)

  #### List reconstruction if only certain levels were selected #####
if(length(datalistselection)!=length(names(datalist))){
##hay que acabar esto la idea es que si no se usan todos los niveles omicos se tome la lista original y se reemplacen en ella los elementos modificados
##determinados por datalistselection
  texts_wizard("\n\nOutput will have the same structure than the original input. Not-processed omic levels will remain unmodified in the output.")
  readline("\nPress [enter] to continue:")
  
  #myprocesseddata <- modifyList(datalist,myprocesseddata,) Do not uncomment, never. 
  datalist[names(myprocesseddata)] <-myprocesseddata
  myprocesseddata <-datalist
  }

#### Data Export ####
class(myprocesseddata)<-"POL"
if(is.null(outputobjectname)==TRUE) {
  outputobjectname <- paste(datalists,"_preprocessed",sep="")
  dog<-as.character(paste("\n\nOutput object name was not set by user. Employing default name: ",outputobjectname, sep=""))
  texts_wizard(dog)
  if(exists(outputobjectname)==T){
  texts_wizard("\nWarning: we found an object with the same name in global environment. Do you want to rename the preprocessed dataset?")
    rename <- utils::select.list(choices = truefalse)
    if(rename == "Yes") rename <- TRUE
    if(rename == "No") rename <- FALSE
   
  if(rename==TRUE){
    texts_wizard("Please provide a new name for your preprocessed dataset and press [enter]")
    outputobjectname <- readline()
    }
  }
}

if(is.character(outputobjectname)!=TRUE) {
  outputobjectname <- paste(datalists,"_preprocessed",sep="")
  dog<-as.character(paste("\n\nOutput object name was wrong. Employing default name: ",paste(datalists,"_preprocessed",sep="")))
  texts_wizard(dog)
  if(exists(outputobjectname)==T){
    texts_wizard("\nWarning: we found an object with the same name in global environment. Do you want to rename the preprocessed dataset?")
    rename <- utils::select.list(choices = truefalse)
    if(rename == "Yes") rename <- TRUE
    if(rename == "No") rename <- FALSE
    
    if(rename==TRUE){
      texts_wizard("Please provide a new name for your preprocessed dataset and press [enter]")
      outputobjectname <- readline()
    }
  }
}
  

if(is.null(treatment2col)==TRUE) treatment2col <- NA
defaultsettings <- c(initialrow=initialrow, initialcolumn=initialcolumn, treatment1col=treatment1col,treatment2col=treatment2col, treatment=treatment)
class(defaultsettings) <- "dsettings"

exportlist<-list(myprocesseddata,defaultsettings)
names(exportlist) <- c(outputobjectname, paste(datalists,"_dsettings",sep=""))
list2env(exportlist,envir = .GlobalEnv)
dog <- as.character(paste("\n",outputobjectname, " has been saved in global environment.", sep=""))
texts_wizard(dog)

texts_wizard("\nDo you want to export the preprocessed dataset to an Excel file?")
answer <- utils::select.list(choices = truefalse)
if(answer == "Yes"){
  texts_wizard("\nPlease provide a name for your XLSX (do not include extension or quotes ;i.e. my_preprocessed_data) file and press [enter]")
  thename<-readline()
  filename <- paste(thename,".xlsx", sep="")
  export_table(myprocesseddata,filename=filename)
  dog <- as.character(paste("\n",filename, " has been saved in your working directory: ", getwd(),sep=""))
  texts_wizard(dog)
}

#### Exporting session log ####
texts_wizard("\npreprocess_Wizard: Duty Accomplished.")
dog <- as.character(paste("Session log, ",sinkfilename, " has been saved in your working directory: ", getwd(),sep=""))
texts_wizard(dog)
sink()


}

checkuserinput <-function(variable,nrow=NULL,ncol=NULL){
  if(is.numeric(variable)==F) stop("FATAL ERROR. Please input an integer number corresponding to desired row/column.")
  if(variable%%1!=0) stop("FATAL ERROR.Please select a correct row/column by introducing an integer. Please re-run.")
  if(is.null(nrow)) if(variable > ncol) stop("FATAL ERROR. The column you have selected is not available in your dataset. Please re-run and select a correct column.")
  if(is.null(ncol)) if(variable > nrow) stop("FATAL ERROR. The row you have selected is not available in your dataset. Please re-run and select a correct row.")
}

