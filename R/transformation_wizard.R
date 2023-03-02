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
#' @name transformation_wizard
#' @title transformation_wizard(): easily transform omic datasets
#' @description A function for quickly and easily transformation of omic datasets
#' @usage transformation_wizard()
#' @details This wizard is aimed to guide you for properly transform your dataset. Please note that this module will only work if you have
#' preprocessed your data with preprocess_wizard() or preprocess_omic_list() functions.
#' @return Transformed data List
#' @author Luis Valledor and Laura Lamelas
#' @seealso \code{\link{importfromexcel}}, \code{\link{preprocess_omic_list}}, \code{\link{preprocess_wizard}}, \code{\link{transformdata}}.
#' @export

transformation_wizard<-function(){
  #### Starting console logging ####
  diahora <-gsub(":","",gsub(" ","_", gsub("-","",Sys.time())))
  diahora <-substr(diahora,3,nchar(diahora)-2)
  sinkfilename <- paste(diahora,"_TransfWzrdLog.txt",sep="")
  sink(sinkfilename, split = T)
  #### Initial text ####
  headings_wizards("Welcome to pRocessomics Data Transformation Wizard")
  texts_wizard("\nThis wizard is aimed to guide you for properly transform your dataset. Please note that this module will only work if you have preprocessed your data with preprocess_wizard() or preprocess_omic_list() functions. Please read 'https://github.com/Valledor/pRocessomics/wiki/Pre-processing-your-data' before continuing.\n")
  texts_wizard("If you have preprocessed the data by yourself please change the class of your datalist to 'POL' before proceeding. But do it at your own risk, since no data structure and quality check will be performed. We strongly recommend using pRocessomics import, annotation, preprocess and transform functions.\n")
  texts_wizard("Press [enter] to continue...")
  invisible(readline())
  
  #### Select datalist ####
  # We ask the user to input one datalist to be used
  datalists <- Filter(function(x) "POL" %in% class(get(x)), ls(envir = .GlobalEnv))
  if(length(datalists)==0) stop("No valid POL objects were found. Please preprocess your data with preprocess_omic_list function or preprocess_wizard.")
  texts_wizard("\nWe have detected the following datasets available. Please enter the number of the dataset do you want to transform and press [enter].")
  texts_wizard("\nIf your dataset is not listed please exit pressing [esc] and re-import your data.")
  datalists<- select.list.edit(datalists, preselect = NULL, multiple = F, title = NULL, graphics = FALSE)
  datalist <-get(datalists,envir=.GlobalEnv)
  dog<-as.character(paste("\nYou have selected ",datalists," dataset. This dataset has the following omic levels: ", paste(names(datalist), collapse=" "),".",sep=""))
  texts_wizard(dog)
  
  
  #### Settings ----
  settings<-selectsettings(originaldatasetname = datalists,datalist = datalist)
  initialrow <- settings[1]
  initialcolumn <- settings[2]
  treatment1col <- settings[3]
  treatment2col <- settings[4]
  treatment <- settings[5]

  
  #### Mini pre-process & treatments -------   
  # Remove all rows before initial
  if(is.null(initialrow)){
    texts_wizard("\n\nWARNING: initial row was not provided. Considering initialrow=1 (Data of first case is in row 1). Please cancel and set proper initialrow if you don't like this setting.\n\n")
    readline(prompt="Press [enter] to continue")
    initialrow=1}
  if(initialrow!=1){
    texts_wizard("\n\nRemoving rows without any relevant data for downstream analyses (rows 1:initialrow).\n")
    dataset <- lapply(dataset, function(x){x <- x[initialrow:nrow(x),]
    return(x)
    })
    initialrow=1
  }
  
  ### Selection of transformation method(s) ####
  texts_wizard("\nAvailable transformations are:")
  tchoices<-c("z","z-centered","Log10","Log10plus1","box-cox-lambda","sqrt","cubic","arcsine","none")
  transf.choose<- select.list.edit(tchoices, preselect = NULL, multiple = T, title = NULL, graphics = FALSE)
  transformation<-transf.choose
  
  #### Datasets to be analyzed ####
  texts_wizard("\nSelect which omic levels you want to process with this settings:")
  datalistselection <- select.list.edit(names(datalist), preselect = NULL, multiple = T, title = NULL, graphics = FALSE)
  finaldatalist <- datalist[datalistselection]
  class(finaldatalist) <- "POL"
  #### Call transformdata ####
  
  texts_wizard("\nLet's shape your data!")
  if(length(transformation)>1){
    dog<-as.character(paste("\nYou have selected ",paste(transformation, collapse=" ")," transformations. Each transformation will be applied independently, generating as many objects as transformation methods you have chosen.",sep=""))
    texts_wizard(dog)
    readline("\nPress [enter] to continue:")
    }
  if(length(datalistselection)!=length(names(datalist))){
    texts_wizard("\n\nYou have decided not transforming all of your omic levels. Outputs will have the same structure than the original input. Not-transformed omic levels will remain unmodified.")
    readline("\nPress [enter] to continue:")  
  }
  for(i in 1:length(transformation))
  {
   
  mytransformeddata<-transformdata(datalist=finaldatalist, initialrow=initialrow, initialcolumn=initialcolumn, 
                                   treatment1col=treatment1col,treatment2col=treatment2col,treatment=treatment,
                                   transf=transf.choose[i])
  
  #### List reconstruction if only certain levels were selected #####
  if(length(datalistselection)!=length(names(datalist))){
   # text have been moved out of the loop
        datalist[names(mytransformeddata)] <-mytransformeddata
  }
  
  #### Saving to environment ####
  outputobjectname <-paste(datalists,"_",transf.choose[i],"_transf", sep="")
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
  class(datalist) <-"POL"
  exportlist<-list(datalist,transf.choose)
  names(exportlist) <- c(outputobjectname, paste(datalists,"_transformsettings",sep=""))
  list2env(exportlist,envir = .GlobalEnv)
  dog <- as.character(paste("\n",outputobjectname, " has been saved in global environment.", sep=""))
  texts_wizard(dog)

  #### Exporting to excel ####
  texts_wizard("\nDo you want to export the transformed dataset to an Excel file?")
  answer <- utils::select.list(choices = truefalse)
  if(answer == "Yes"){
    texts_wizard("\nPlease provide a name for your XLSX (do not include extension;i.e. .xlsx) file and press [enter]")
    thename<-readline()
    filename <- paste(thename,".xlsx", sep="")
    export_table(datalist,filename=filename)
    dog <- as.character(paste("\n",filename, " has been saved in your working directory: ", getwd(),sep=""))
    texts_wizard(dog)
  }
  }

  #### Exporting session log ####
  texts_wizard("\nTransformation_Wizard: Duty Accomplished.")
  dog <- as.character(paste("Session log, ",sinkfilename, " has been saved in your working directory: ", getwd(),sep=""))
  texts_wizard(dog)
  sink()
  
}