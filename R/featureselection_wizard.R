#
# Copyright 2019 Luis Valledor / Laura Lamelas
# Last revision 04.09.2019
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
#' @name featureselection_wizard
#' @title featureselection_wizard() easily filter more significant variablest
#' @description This function is a wizard to interactively help the user to perform the feature selection on preprocessed data
#' @usage featureselection_wizard()
#' @seealso featureselect
#' @details The objective of \code{featureselection_wizard} is removing those variables who do (not) 
#' greatly change between treatments, 
#' known to later introduce noise in multivariate analyses reducing the power of these approaches. 
#' This function requires a \code{POL} object, previously defined by using \code{preprocess_omic_list} and, 
#' if desired, \code{transformdata}. Two different approaches are implemented within this 
#' function: first one is based on selecting those variables that have inter-quartilic ranges, 
#' IQR, greater than average within all studied variables. In this case the \code{threshold} i
#' ndicates the \% above average IQR and only variables with greater or equal IQR will be selected. 
#' The second strategy is based on applying ANOVA or Kruskal Wallis and selected those variables 
#' exhibiting a p or q value smaller or equal than \code{threshold}.
#' @return returns a "POL" object with the filtered data
#' @author Luis Valledor and Laura Lamelas
#' @seealso featureselection
#' @export
#'


featureselection_wizard<-function(){
  
  #### Starting console logging ####
  diahora <-gsub(":","",gsub(" ","_", gsub("-","",Sys.time())))
  diahora <-substr(diahora,3,nchar(diahora)-2)
  sinkfilename <- paste(diahora,"_FeatSelLog.txt",sep="")
  sink(sinkfilename, split = T)
  
  #### Initial texts ####
  headings_wizards("Welcome to pRocessomics Feature Selection Wizard")
  texts_wizard("\nThis wizard is aimed to guide you for properly filtering noisy variables in your dataset. Please note that this module will only work if you have preprocessed your data with preprocess_wizard or preprocess_omic_list functions. Please read 'https://github.com/Valledor/pRocessomics/wiki/Pre-processing-your-data' before continuing.")
  texts_wizard("\nIf you have preprocessed the data by yourself please change the class of your datalist to 'POL' before proceeding. But do it at your own risk, since no data structure and quality check will be performed. We strongly recommend using pRocessomics import, annotation, preprocess and transform functions.")
  texts_wizard("\nPress [enter] to continue...")
  invisible(readline())
  
  
  #### Select datalist ####
  # We ask the user to input one datalist to be used
  datalists <- Filter(function(x) "POL" %in% class(get(x)), ls(envir = .GlobalEnv))
  if(length(datalists)==0) stop("No valid POL objects were found. Please preprocess your data with preprocess_omic_list function or preprocess_wizard.")
  texts_wizard("\nWe have detected the following datasets available. Please enter the number of the dataset do you want to filter and press [enter].")
  texts_wizard("\nIf your dataset is not listed please exit pressing [esc] and re-import your data.")
  datalists<- select.list.edit(datalists, preselect = NULL, multiple = F, title = NULL, graphics = FALSE)
  datalist <-get(datalists,envir=.GlobalEnv)
  dog<-as.character(paste("\nYou have selected ",datalists," dataset. This dataset has the following omic levels: ", paste(names(datalist), collapse=" "),".",sep=""))
  texts_wizard(dog)
  
  
  
    #### Settings ----
  #### Settings ----
  settings<-selectsettings(originaldatasetname = datalists,datalist = datalist)
  initialrow <- settings[1]
  initialcolumn <- settings[2]
  treatment1col <- settings[3]
  treatment2col <- settings[4]
  treatment <- settings[5]
  parametric = T 
  stat = "p"
  abovebelow = "Below"
  
  #### Mini pre-process & treatments -------   
  # Remove all rows before initial
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
  
  ### Selection of feature selection method(s) ####
  texts_wizard("\nThere are three available methods for droping noisy or constant variables of your dataset: InterQuartile Range (IQR), variables with an IQR Above/Below threshold will be removed; p-value (Stat), variables showing a p-value (ANOVA or Kruskal-Wallis) greater/lower than threshold will be removed; or coefficient of variation (VarCoef), variables showing a coefficient of variation greater/lower than threshold will be removed. Please check 'https://github.com/Valledor/pRocessomics/wiki/Pre-processing-your-data' for further information. \nPlease select which filtering method(s) you want to apply")
  tchoices<-c("IQR","Stat","VarCoef")
  featuresel.choose<- select.list.edit(tchoices, preselect = NULL, multiple = T, title = NULL, graphics = FALSE)
  featureselected<-featuresel.choose
  
  #### Datasets to be analyzed ####
  texts_wizard("\nSelect which omic levels you want to process with this settings:")
  datalistselection <- select.list.edit(names(datalist), preselect = NULL, multiple = T, title = NULL, graphics = FALSE)
  finaldatalist <- datalist[datalistselection]
  class(finaldatalist) <- "POL"
  
  #### Call featureselect ####
  texts_wizard("\nLet's filter your data!")
  
  if(length(featureselection)>1){
    dog<-as.character(paste("\nYou have selected ",paste(featureselected, collapse=" ")," methods. Each filtering method will be applied independently, generating as many objects as methods you have chosen.",sep=""))
    texts_wizard(dog)
    readline("\nPress [enter] to continue:")
  }
  if(length(datalistselection)!=length(names(datalist))){
    texts_wizard("\n\nYou have decided not filtering all of your omic levels. Unfiltered outputs will remain unmodified.")
    readline("\nPress [enter] to continue:")  
  }
  for(i in 1:length(featureselected)) {
  if(featureselected[i]=="IQR"){
    dog<-"\nYou have selected IQR. First we will calculate the average IQR of your dataset, and then and multiply it by your threshold. All variables with an IQR above/below (you will also choose that) this value will be kept."
    texts_wizard(dog)
    dog<-"Please provide threshold and press [enter]: \n"
    texts_wizard(dog)
    threshold <- readline()
    threshold<-as.numeric(threshold)
    if(threshold<=-100|threshold>=100) stop("Threshold out of limits. Please enter a number between 1 and 99")
    threshold=threshold/100
    dog<-"\nDo you want to keep variables above or below this threshold?"
    texts_wizard(dog)
    tchoices<-c("Above","Below")
    abovebelow<- select.list.edit(tchoices, preselect = NULL, multiple = F, title = NULL, graphics = FALSE)
  }
  
    if(featureselected[i]=="Stat"){
      dog<-"\nYou have selected Stat. First you will choose p-value threshold, then if you want to get these p-values employing a parametric (ANOVA) or non-parametric (Kruskal-Wallis) approach, and finally if you want to keep those variables above/below threshold."
      texts_wizard(dog)
      dog<-"Please provide threshold for p-values [0.001-0.999] and press [enter]: \n"
      texts_wizard(dog)
      threshold <- readline()
      threshold<-as.numeric(threshold)
      if(threshold<=0.0009|threshold>=1) stop("Threshold out of limits. Please enter a number between 1 and 99")

      texts_wizard("\n Do you want to employ parametric (ANOVA) analysis?")
      parachoice <- c(TRUE,FALSE)
      parametric<-select.list.edit(parachoice, preselect = NULL, multiple = F, title = NULL, graphics = FALSE)
      statchoice<-c("p","q")
      texts_wizard("\n Which stat would you like to use, p or q value?")
      stat<-select.list.edit(statchoice, preselect = NULL, multiple = F, title = NULL, graphics = FALSE)

      texts_wizard("\nDo you want to keep variables above or below this threshold?")
      tchoices<-c("Above","Below")
      abovebelow<- select.list.edit(tchoices, preselect = NULL, multiple = F, title = NULL, graphics = FALSE)
    }  
  
    if(featureselected[i]=="VarCoef"){
      dog<-"\nYou have selected coefficient of variation (CV). First we will calculate the average CV, and then we will filter those variables which CV is xx% greater/lower thant the average. First you will choose threshold XX%, and then if you want to keep those variables above/below threshold."
      texts_wizard(dog)
      dog<-"Please provide threshold [1-99] and press [enter]: \n"
      texts_wizard(dog)
      threshold <- readline()
      threshold<-as.numeric(threshold)
      if(threshold<=0|threshold>=100) stop("Threshold out of limits. Please enter a number between 1 and 99")
      threshold=threshold/100
      dog<-"\nDo you want to keep variables above or below this threshold?"
      texts_wizard(dog)
      tchoices<-c("Above","Below")
      abovebelow<- select.list.edit(tchoices, preselect = NULL, multiple = F, title = NULL, graphics = FALSE)
    }
        
  filtereddata<-featureselection(datalist=finaldatalist, initialrow=initialrow, initialcolumn=initialcolumn, 
                                      treatment1col=treatment1col,treatment2col=treatment2col,treatment=treatment,
                                      method = featureselected[i], threshold = threshold, parametric=parametric, stat = stat,abovebelow=abovebelow)
  
  #### List reconstruction if only certain levels were selected #####
  if(length(datalistselection)!=length(names(datalist))){
    # text have been moved out of the loop
    datalist[names(filtereddata)] <-filtereddata
  } else{
  datalist <- filtereddata
  }
  #### Saving to environment ####
  outputobjectname <-paste(datalists,"_",featureselected[i],"_featsel", sep="")
  if(exists(outputobjectname)==T){
    texts_wizard(paste("\nWarning: we found an object with", outputobjectname, "name in global environment. Do you want to rename the preprocessed dataset?",sep=" "))
    rename <- utils::select.list(choices = truefalse)
    if(rename == "Yes") rename <- TRUE
    if(rename == "No") rename <- FALSE
    
    if(rename==TRUE){
      texts_wizard("Please provide a new name for your preprocessed dataset and press [enter]")
      outputobjectname <- readline()
    }
  }
  class(datalist) <-"POL"
  exportlist<-list(datalist,featuresel.choose)
  names(exportlist) <- c(outputobjectname, paste(datalists,"_featselsettings",sep=""))
  list2env(exportlist,envir = .GlobalEnv)
  dog <- as.character(paste("\n",outputobjectname, " has been saved in global environment.", sep=""))
  texts_wizard(dog)
  
  #### Exporting to excel ####
  texts_wizard("\nDo you want to export the filtered dataset to an Excel file?")
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
  texts_wizard("\nFeatureSelection_Wizard: Duty Accomplished.")
  dog <- as.character(paste("Session log, ",sinkfilename, " has been saved in your working directory: ", getwd(),sep=""))
  texts_wizard(dog)
  sink()
}