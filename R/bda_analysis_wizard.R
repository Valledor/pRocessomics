#
# Copyright 2020 Luis Valledor / Laura Lamelas
# Last revision 05.05.2020
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
#' @name bda_analysis_wizard
#' @title block-DA wizard
#' @description A function for easily integrate multi-omics data (block discriminant analysis)
#' @usage bda_analysis_wizard("outputobjectname")
#' @param outputobjectname Character, the quoted name of the result of the analysis
#' @author Luis Valledor and Laura Lamelas
#' @return Block DA results, a "bdaanalysis" class object
#' @importFrom stats prcomp
#' @importFrom utils select.list
#' @export

bda_analysis_wizard <- function(outputobjectname=NULL){
  #### Starting console logging ####
  diahora <-gsub(":","",gsub(" ","_", gsub("-","",Sys.time())))
  diahora <-substr(diahora,3,nchar(diahora)-2)
  sinkfilename <- paste(diahora,"_BDAWzrdLog.txt",sep="")
  sink(sinkfilename, split = T)
  
  #### Initial texts ####
  headings_wizards("Welcome to pRocessomics block-DA Wizard")
  texts_wizard("\n\nThis wizard is aimed to guide you for properly perform a block-DA on your datasets, including model tuning and performance estimation. Please note that this module will only work if data has been already preprocessed with preprocess_wizard() or preprocess_omic_list() functions. Please read 'https://github.com/Valledor/pRocessomics/wiki/Multivariate' before continuing.\n")
  texts_wizard("\nIf you have preprocessed the data by yourself please change the class of your datalist to 'POL' before proceeding. But do it at your own risk, no data structure and quality check will be performed. We strongly recommend using pRocessomics import, annotation, preprocess and transform functions.")
  texts_wizard("\nMost of the function employed by this wizard are part of MixOmics package, please visit 'http://mixomics.org/' for more information and please do not forget citing the work of our colleagues.")
  texts_wizard("\nPress [enter] to continue...")
  invisible(readline())
  
  #### Select datalist ####
  datalists <- Filter(function(x) "POL" %in% class(get(x)), ls(envir = .GlobalEnv))
  if(length(datalists)==0) stop("No valid POL objects were found. Please preprocess your data with preprocess_omic_list function or preprocess_wizard.")
  texts_wizard("\nWe have detected the following datasets available. Please enter the number of the dataset do you want to choose to perform the block-DA analysis and press [enter].")
  texts_wizard("\nIf your dataset is not listed please exit pressing [esc]  and re-import your data.")
  
  datalists <- select.list.edit(datalists, preselect = NULL, multiple = F, title = NULL, graphics = FALSE)
  writeLines(datalists, con = stdout(), sep = "\n", useBytes = FALSE) 
  
  datalist <- get(datalists,envir=.GlobalEnv)
  datasetnames <- names(datalist)
  dog <- as.character(paste("You have selected",datalists,"dataset. This dataset has the following omic levels:\n"))
  texts_wizard(dog)
  texts_wizard(paste(datasetnames,collapse=", "))
  if(length(datasetnames)<2) stop("FATAL ERROR: Your dataset must have, at least, two diferent omic levels")
  tempdatasetname <- datalists
  
  #### Select omic levels----
  datasetnames<-names(datalist)
  texts_wizard("\nWe have detected the following levels in your data. Please select at least two levels to be used and press [enter].\n")
  datasetnames <- c(datasetnames)
  omiclevels <- utils::select.list(datasetnames, preselect = NULL, multiple = T, title = NULL, graphics = FALSE)
  writeLines(omiclevels, con = stdout(), sep = "\n", useBytes = FALSE) 
  if(length(omiclevels)<2) stop("FATAL ERROR: You must choose at least two levels")

  
  #### dSettings and treatment selection ----
  settings<-selectsettings(originaldatasetname = datalists,datalist = datalist)
  initialrow <- settings[1]
  initialcolumn <- settings[2]
  treatment1col <- settings[3]
  treatment2col <- settings[4]
  
  # Treatment definition 
  if(!is.na(treatment2col)){
    texts_wizard("\nPlease select which treatment you want to use as factor:\n")
    treatment<-select.list.edit(c("Treatment 1", "Treatment 2","Combination of Treatments 1 and 2"), preselect = NULL, multiple = F, title = NULL, graphics = FALSE)
    writeLines(treatment, con = stdout(), sep = "\n", useBytes = FALSE) 
    
    if(treatment=="Treatment 1") treatment<-1
    if(treatment=="Treatment 2") treatment<-2
    if(treatment=="Combination of Treatments 1 and 2") treatment <-3
  }
  if(is.na(treatment2col)){
    treatment=1
    treatment2col <-treatment1col}
  
  #### Select and add annotations ----
  texts_wizard("\nDo you want to add annotations to this analysis?")
  an <- utils::select.list(choices = truefalse)
  if(an == "Yes") an <- TRUE
  if(an == "No") an <- FALSE
  writeLines(an, con = stdout(), sep = "\n", useBytes = FALSE) 
  
  if(an==T){
    annotation <- selectannotation()
    class(annotation)<-"data.frame"
    annotation<-as.data.frame(annotation)
  }
  else{
    annotation<-NULL
  }
  
  ####Autotune, performance, and bda_analysis ----
  texts_wizard("\nBDA model allows to define the number of variables to be kept after analysis. Depending on the number of retained variables the model will have different performances (see 'https://github.com/Valledor/pRocessomics/wiki/Multivariate'). Do you want to autotune the model (yes) or employ all variables (no)?")
  autotune <- utils::select.list(choices = truefalse)
  writeLines(autotune, con = stdout(), sep = "\n", useBytes = FALSE) 
  if(autotune=="Yes") autotune<-T
  if(autotune=="No") autotune<-F
 
  bda_results <- bda_analysis(datalist, omiclevel=omiclevels, annotation, initialrow, initialcolumn, treatment1col, treatment2col, treatment, autotune =autotune)
  texts_wizard("Press [enter] to continue")
  invisible(readline())
  
  #### Naming outputs ####
  if(is.null(outputobjectname)==TRUE) {
    outputobjectname <- paste(datalists,"_bda",sep="")
    dog<-as.character(paste("\n\nOutput object name was not set by user. Employing default name: ",paste(datalists,"_bda",sep="")))
    texts_wizard(dog)
    if(exists(outputobjectname)==T){
      texts_wizard("\nNotice that there is an object in global environment with the same name and **it will be replaced**")
      texts_wizard("\nDo you want to change this name and provide a unique identifier?")
      rename <- utils::select.list(choices = truefalse)
      if(rename == "Yes") rename <- TRUE
      if(rename == "No") rename <- FALSE
      writeLines(rename, con = stdout(), sep = "\n", useBytes = FALSE) 
      
      if(rename==TRUE){
        newsuffix<-readline("Please provide a suffix for your dataset _")
        outputobjectname <- paste(datalists,newsuffix,sep="_")
        writeLines(outputobjectname, con = stdout(), sep = "\n", useBytes = FALSE) 
        
      }
    }
  }
  
  if(is.character(outputobjectname)!=TRUE) {
    outputobjectname <- paste(datalists,"_bda",sep="")
    dog<-as.character(paste("\n\nOutput object name was wrong. Employing default name: ",paste(tempdatasetname,"_bda",sep="")))
    texts_wizard(dog)
    if(exists(outputobjectname)==T){
      texts_wizard("\nnotice that there is an object in global environment with the same name and **it will be replaced**")
      texts_wizard("\nDo you want to change this name and provide a unique identifier?")
      rename <- utils::select.list(choices = truefalse)
      if(rename == "Yes") rename <- TRUE
      writeLines(rename, con = stdout(), sep = "\n", useBytes = FALSE) 
      
      if(rename==TRUE){
        newsuffix<-readline("Please provide a suffix for your dataset _")
        outputobjectname <- paste(datalists,newsuffix,sep="_")
        writeLines(newsuffix, con = stdout(), sep = "\n", useBytes = FALSE) 
        
      }
    }
  }
  
  
  #### Export the result to R and to excel----
  exportlist <- list(bda_results)
  names(exportlist) <- c(outputobjectname)
  list2env(exportlist,envir = .GlobalEnv)
  dog <- as.character(paste(outputobjectname,"has been exported successfuly to your global environment",sep=" "))
  texts_wizard(dog)
  
  texts_wizard("\nDo you want to save block-DA results in .xlsx?")
  answer <- utils::select.list(choices = truefalse)
  writeLines(answer, con = stdout(), sep = "\n", useBytes = FALSE) 
  
  if(answer == "Yes"){
    texts_wizard("Saving...")
    diahora <-gsub(":","",gsub(" ","_", gsub("-","",Sys.time())))
    diahora <-substr(diahora,3,nchar(diahora)-2)
    tablefilename<-paste(diahora,"_bda_results.xlsx",sep="")
    export_table(result_list = bda_results, filename = tablefilename)
    texts_wizard(paste("\nResults saved in your working directory,",paste(getwd(),", as",sep=""), tablefilename,".", sep=" "))
  }
  
  
  #### Guided plotting ####
  texts_wizard("\nAre your ready to plot block-DAS results?")
  
  answer <- utils::select.list(choices = truefalse)
  writeLines(answer, con = stdout(), sep = "\n", useBytes = FALSE) 
  
  while (answer == "Yes"){
    texts_wizard("\nPlease select the results you want to plot")
    texts_wizard("Distance, represents samples according to their scores to selected components")
    texts_wizard("Distance_Ellipse, Distance plot with confidence ellipses according to sample groups")
    texts_wizard("Topscoring, top scoring variables in the selected component")
    texts_wizard("Cim, Clustered Image Map summary of the analysis")
    texts_wizard("Network, representation of the relation between variables and selection vector")
    texts_wizard("Circos, representation of vriable correlation through omic levels")
    #plottype
    bdaplotstype<-c("Distance","Distance_Ellipse", "Topscoring", "Cim", "Network","Circos")
    plottype <- utils::select.list(choices = bdaplotstype)
    writeLines(plottype, con = stdout(), sep = "\n", useBytes = FALSE) 
    
    #Components
    CompX=1
    CompY=2
    cutoff=0
    confidence= 0.9
    treatment = 1
    ### treatment ####
    if (plottype %in% c("Distance", "Distance_Ellipse", "Cim")){
      texts_wizard("\nPlease select which treatment you want to use as factor:\n")
      treatment<-select.list.edit(c("Treatment 1", "Treatment 2","Combination of Treatments 1 and 2"), preselect = NULL, multiple = F, title = NULL, graphics = FALSE)
      writeLines(treatment, con = stdout(), sep = "\n", useBytes = FALSE) 
      
      if(treatment=="Treatment 1") treatment<-1
      if(treatment=="Treatment 2") treatment<-2
      if(treatment=="Combination of Treatments 1 and 2") treatment <-3
    }
    #Ellipse#####
    if(plottype == "Distance_Ellipse"){
      texts_wizard("Please type the confidence value for the ellipses. (type 0.95 for 95% confidence ellipses)")
      confidence <- as.numeric(readline())
      writeLines(confidence, con = stdout(), sep = "\n", useBytes = FALSE) 
      
    }
    #Annotations####
    zzz<-F
    if(plottype %in% c("Topscoring","Network")){
      anotlists <- Filter(function(x) "pRoAnnot" %in% class(get(x)), ls(envir = .GlobalEnv))
      if(an==T) {
        texts_wizard("Do you want to use annotations with your plot?")
        zzz<-utils::select.list(truefalse)
        writeLines(zzz, con = stdout(), sep = "\n", useBytes = FALSE) 
        
        if(zzz=="Yes") zzz<-T
        if(zzz=="No") zzz<-F
        if(zzz==TRUE) {
          if(plottype=="Composite") plottype<-"Composite_Annotated"
        }
      }
    }
    #Fortopscoring ####
    fortopscoring <- c(1,15,"abs")
    if(plottype =="Topscoring"){
      texts_wizard("Which component number do you want to analyze?")
      a<-(readline())
      writeLines(a, con = stdout(), sep = "\n", useBytes = FALSE) 
      
      texts_wizard("How many of top scoring variables do you want to plot?")
      b<-(readline())
      writeLines(b, con = stdout(), sep = "\n", useBytes = FALSE) 
      
      texts_wizard("There are currently two modes for screening the variables to display: abs for absolute value scores and posneg for positive and negative loading analyzed separately, in which case the top x positive and negative variables will be displayed.")
      texts_wizard("Which mode do you want to use?")
      opts <- c("abs","posneg")
      c <- utils::select.list(opts)
      writeLines(c, con = stdout(), sep = "\n", useBytes = FALSE) 
      
      fortopscoring <-c(a,b,c)
      
    } 
    #cutoff####
    if(plottype %in% c("Network", "Circos")){
      texts_wizard("Please select the correlation value cutoff for edge consideration. ")
      cutoff <- readline("\nIntroduce a value between 0 and 1 and press [enter]:")
      writeLines(cutoff, con = stdout(), sep = "\n", useBytes = FALSE) 
      
      if(cutoff!="") {
        cutoff <- as.numeric(cutoff)
        if((cutoff>=0&&cutoff<=1)==FALSE){
          texts_wizard("FATAL ERROR: cutoff should be a decimal value between 0 and 1")
          next()
        } 
      }  
    } 
   
    #### Call bda_plot ####
    plot.object <- bda_plot(bda_list = bda_results,treatment = treatment,plottype = plottype,fortopscoring = fortopscoring, useannot = zzz,cutoff = cutoff, confidence = confidence)
    print(plot.object)
    
    #### Exporting bda_plot object ####  
    texts_wizard("Do you want to export the current plot?")
    ansplot <- utils::select.list(choices = truefalse)
    writeLines(ansplot, con = stdout(), sep = "\n", useBytes = FALSE) 
    
    if(ansplot=="Yes"){
      texts_wizard("Saving...")
      diahora <-gsub(":","",gsub(" ","_", gsub("-","",Sys.time())))
      diahora <-substr(diahora,3,nchar(diahora)-2)
      plotfilename<-paste(diahora,"_bda_",plottype,"_plot.pdf",sep="")
      export_plot(plot.object, plotfilename)
      
      texts_wizard(paste("\nPlot saved in your working directory,",paste(getwd(),", as",sep=""),plotfilename, sep=" "))
    }
    texts_wizard("\nDo you want to try any other plot?")
    answer <- utils::select.list(choices = truefalse)
    writeLines(answer, con = stdout(), sep = "\n", useBytes = FALSE) 
    
  }
  
  #### Exporting session log ####
  texts_wizard("\nbda_analysis_Wizard: Duty Accomplished.")
  dog <- as.character(paste("Session log, ",sinkfilename, " has been saved in your working directory: ", getwd(),sep=""))
  texts_wizard(dog)
  sink()
  #sink()
}




