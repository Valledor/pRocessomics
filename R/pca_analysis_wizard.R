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
#' @name pca_analysis_wizard
#' @title PCA wizard
#' @description A function for easily perform a principal component analysis
#' @usage pca_analysis_wizard("outputobjectname")
#' @param outputobjectname Character, the quoted name of the result of the analysis
#' @author Luis Valledor and Laura Lamelas
#' @return The PCA results, a "pcaanalysis" class object
#' @importFrom stats prcomp
#' @importFrom utils select.list
#' @export

pca_analysis_wizard <- function(outputobjectname=NULL){
  #### Starting console logging ####
  diahora <-gsub(":","",gsub(" ","_", gsub("-","",Sys.time())))
  diahora <-substr(diahora,3,nchar(diahora)-2)
  sinkfilename <- paste(diahora,"_PCAWzrdLog.txt",sep="")
  sink(sinkfilename, split = T)
  
  
  #### Initial texts ####
  headings_wizards("Welcome to pRocessomics PCA Wizard")
  texts_wizard("\n\nThis wizard is aimed to guide you for properly perform a PCA on your datasets. Please note that this module will only work if data has been already preprocessed with preprocess_wizard() or preprocess_omic_list() functions.Please read 'https://github.com/Valledor/pRocessomics/wiki/Multivariate' before continuing.\n")
  texts_wizard("\nIf you have preprocessed the data by yourself please change the class of your datalist to 'POL' before proceeding. But do it at your own risk, no data structure and quality check will be performed. We strongly recommend using pRocessomics import, annotation, preprocess and transform functions.")
  texts_wizard("\nPress [enter] to continue...")
  invisible(readline())
  
  #### Select datalist ####
  datalists <- Filter(function(x) "POL" %in% class(get(x)), ls(envir = .GlobalEnv))
  if(length(datalists)==0) stop("No valid POL objects were found. Please preprocess your data with preprocess_omic_list function or preprocess_wizard.")
  texts_wizard("\nWe have detected the following datasets available. Please enter the number of the dataset do you want to choose to perform the PCA analysis and press [enter].")
  texts_wizard("\nIf your dataset is not listed please exit pressing [esc]  and re-import your data.")
  
  datalists<- select.list.edit(datalists, preselect = NULL, multiple = F, title = NULL, graphics = FALSE)
  datalist <-get(datalists,envir=.GlobalEnv)
  dog<-as.character(paste("You have selected",datalists,"dataset. This dataset has the following omic levels:\n"))
  texts_wizard(dog)
  texts_wizard(paste((names(datalist)),collapse=", "))
  tempdatasetname <- datalists
  
  #### Select omic levels----
  datasetnames<-names(datalist)
  if(length(datasetnames)>1){
    texts_wizard("\nWe have detected the following levels in your data. Please enter the number(s) of the level(s) you want to process and press [enter]\n")
    datasetnames <- c(datasetnames)
    omiclevel <- utils::select.list(datasetnames, preselect = NULL, multiple = T, title = NULL, graphics = FALSE)
  } else {
    omiclevel <- datasetnames
  }
  
  #### dSettings and treatment selection ----
  settings<-selectsettings(originaldatasetname = datalists,datalist = datalist)
  initialrow <- settings[1]
  initialcolumn <- settings[2]
  treatment1col <- settings[3]
  treatment2col <- settings[4]
  
  # Treatment definition 
  if(!is.na(treatment2col)){
    texts_wizard("\nPlease select which treatment you want to use for splitting your data (needed for later scoreplots).\n")
    treatment<-select.list.edit(c("Treatment 1", "Treatment 2","Combination of Treatments 1 and 2"), preselect = NULL, multiple = F, title = NULL, graphics = FALSE)
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
  
  if(an==T){
    annotation<-selectannotation()
    class(annotation)<-"data.frame"
    annotation<-as.data.frame(annotation)
    
  }
  else{
    annotation<-NULL
  }
  
  pca_results <- pca_analysis(datalist, annotation, initialrow, initialcolumn, treatment1col, treatment2col, treatment,omiclevel)
  texts_wizard("Press [enter] to continue")
  invisible(readline())
  
  #### Naming outputs ####
  if(is.null(outputobjectname)==TRUE) {
    outputobjectname <- paste(datalists,"_pca",sep="")
    dog<-as.character(paste("\n\nOutput object name was not set by user. Employing default name: ",paste(datalists,"_pca",sep="")))
    texts_wizard(dog)
    if(exists(outputobjectname)==T){
      texts_wizard("\nNotice that there is an object in global environment with the same name and **it will be replaced**")
      texts_wizard("\nDo you want to change this name and provide a unique identifier?")
      rename <- utils::select.list(choices = truefalse)
      if(rename == "Yes") rename <- TRUE
      if(rename == "No") rename <- FALSE
      
      if(rename==TRUE){
        newsuffix<-readline("Please provide a suffix for your dataset _")
        outputobjectname <- paste(datalists,newsuffix,sep="_")
      }
    }
  }
  
  if(is.character(outputobjectname)!=TRUE) {
    outputobjectname <- paste(datalists,"_pca",sep="")
    dog<-as.character(paste("\n\nOutput object name was wrong. Employing default name: ",paste(tempdatasetname,"_pca",sep="")))
    texts_wizard(dog)
    if(exists(outputobjectname)==T){
      texts_wizard("\nnotice that there is an object in global environment with the same name and **it will be replaced**")
      texts_wizard("\nDo you want to change this name and provide a unique identifier?")
      rename <- utils::select.list(choices = truefalse)
      if(rename == "Yes") rename <- TRUE
      
      if(rename==TRUE){
        newsuffix<-readline("Please provide a suffix for your dataset _")
        outputobjectname <- paste(datalists,newsuffix,sep="_")
      }
    }
  }
  
  
  #### Export the result to R and to excel----
  exportlist <- list(pca_results)
  names(exportlist) <- c(outputobjectname)
  list2env(exportlist,envir = .GlobalEnv)
  dog <- as.character(paste(outputobjectname,"has been exported successfuly to your global environment",sep=" "))
  texts_wizard(dog)
  
  texts_wizard("\nDo you want to save PCA results in .XLSX?")
  answer <- utils::select.list(choices = truefalse)
  if(answer == "Yes"){
    texts_wizard("Saving...")
    diahora <-gsub(":","",gsub(" ","_", gsub("-","",Sys.time())))
    diahora <-substr(diahora,3,nchar(diahora)-2)
    tablefilename<-paste(diahora,"_pca_results.xlsx",sep="")
    export_table(result_list = pca_results, filename = tablefilename)
    texts_wizard(paste("\nResults saved in your working directory,",paste(getwd(),", as",sep=""), tablefilename,".", sep=" "))
  }
  
  
  #### Guided plotting ####
  texts_wizard("\nAre your ready to plot PCA results?")
  
  answer <- utils::select.list(choices = truefalse)
  while (answer == "Yes"){
    texts_wizard("\nPlease select the results you want to plot")
    texts_wizard("Scree, showing the variance explained by each component, and also the accumulative.")
    texts_wizard("Score, represents samples according to their scores to selected PCs")
    texts_wizard("Score_Ellipse, a Score plot with confidence ellipses by sample groups")
    texts_wizard("Topscoring, top scoring variables in the selected principal component")
    texts_wizard("Var, represents variables according to their scores to selected PCs")
    texts_wizard("Biplot, PCA biplot showing samples and variables, Score + Var")
    #texts_wizard("Composite: combination of the most representative plots (Score, Var and Biplot)")
    #plottype
    pcaplotstype<-c("Scree", "Score","Score_Ellipse", "Topscoring",  "Var","Biplot", "Composite")
    plottype <- utils::select.list(choices = pcaplotstype)
    #treatment
    treatment = 1
    if(plottype %in% c("Score","Score_Ellipse","Topscoring","Biplot")){
      if(treatment2col != treatment1col){
        texts_wizard("\nPlease select which treatment you want to use for splitting your data.\n")
        treatment<-select.list(c("Treatment 1", "Treatment 2","Combination of Treatments 1 and 2","Treatment 1 colors, Treatment 2 symbols"), preselect = NULL, multiple = F, title = NULL, graphics = FALSE)
        if(treatment=="Treatment 1") treatment <- 1
        if(treatment=="Treatment 2") treatment <- 2
        if(treatment=="Combination of Treatments 1 and 2") treatment <- 3
        if(treatment=="Treatment 1 for colors and Treatment 2 for symbols") treatment <- 4
      }
      if(treatment2col==treatment1col) treatment <- 1
      
    }
    #Components
    ncomp=length(pca_results$pca$sdev)
    CompX=1
    CompY=2
    CompZ=NULL
    if(plottype %in% c("Score","Score_Ellipse","Biplot","Var")){
      texts_wizard("Now you will select which components will be plot")
      texts_wizard("Component to be plotted in X axis (type 1 for PC1, 2 for PC2, ...)")
      CompX <- as.numeric(readline())
      while(CompX > ncomp) {
        dog <-as.character(paste("Watch out!, There are only ", ncomp, "in the current PCA result object. Please type again the component number to be plotted in X axis"))
        texts_wizard(dog)
        CompX <- as.numeric(readline())
      }
      texts_wizard("Component to be plotted in Y axis (type 1 for PC1, 2 for PC2, ...)")
      CompY <- as.numeric(readline())
      while(CompY == CompX) {
        dog <-as.character(paste("You have already chosen this component. Please select a different component and press [enter]."))
        texts_wizard(dog)
        CompY <- as.numeric(readline())
      }
      while(CompY > ncomp) {
        dog <-as.character(paste("Watch out!, There are only ", ncomp, "in the current PCA result object. Please type again the component number to be plotted in Y axis"))
        texts_wizard(dog)
        CompY <- as.numeric(readline())
      }
    }
    if(plottype == "Score"){
      texts_wizard("3D plot?")
      threeD<-utils::select.list(truefalse)
      if(threeD=="Yes") {
        texts_wizard("Component to be plotted in Z axis (type 1 for PC1, 2 for PC2, ...)")
        CompZ <- as.numeric(readline())
        while(CompZ == CompX|CompZ==CompY) {
          dog <-as.character(paste("You have already chosen this component. Please select a different component and press [enter]."))
          texts_wizard(dog)
          CompZ <- as.numeric(readline())
        }
        while(CompZ > ncomp) {
          dog <-as.character(paste("Watch out!, There are only ", ncomp, "in the current PCA result object. Please type again the component number to be plotted in Z axis"))
          texts_wizard(dog)
          CompZ <- as.numeric(readline())
        }
      }
      if(threeD=="No") CompZ <- NULL
    }
    confidence = 0.9
    if(plottype == "Score_Ellipse"){
      texts_wizard("Please type the confidence value for the ellipses. (type 0.95 for 95% confidence ellipses)")
      confidence <- as.numeric(readline())
    }
    #Annotations
    zzz<-F
    if(plottype %in% c("Topscoring","Biplot","Var")){
      anotlists <- Filter(function(x) "pRoAnnot" %in% class(get(x)), ls(envir = .GlobalEnv))
      if(an==T) {
        texts_wizard("Do you want to use annotations with your plot?")
        zzz<-utils::select.list(truefalse)
        if(zzz=="Yes") zzz<-T
        if(zzz=="No") zzz<-F
      }
    }
    #Fortopscoring
    fortopscoring <- c(1,15,"abs")
    if(plottype =="Topscoring"){
      texts_wizard("Which component number do you want to analyze?")
      a<-(readline())
      texts_wizard("How many of top scoring variables do you want to plot?")
      b<-(readline())
      texts_wizard("Which mode do you want to use?")
      texts_wizard("the mode to compute abs for absolute value and posneg for positive and negative loading analyzed separately")
      opts <- c("abs","posneg")
      c <- utils::select.list(opts)
      fortopscoring <-c(a,b,c)
    }  
    
    #### Call pca_plot ####
    plot.object<-pca_plot(pca_list = pca_results,treatment = treatment, compX = CompX, compY = CompY,compZ = CompZ,plottype = plottype,fortopscoring = fortopscoring, useannot = zzz, confidence= confidence)
    print(plot.object)
    
    #### Exporting pca_plot object ####  
    texts_wizard("Do you want to export the current plot?")
    ansplot <- utils::select.list(choices = truefalse)
    if(ansplot=="Yes"){
      texts_wizard("Saving...")
      diahora <-gsub(":","",gsub(" ","_", gsub("-","",Sys.time())))
      diahora <-substr(diahora,3,nchar(diahora)-2)
      plotfilename<-paste(diahora,"_pca_",plottype,"_plot.pdf",sep="")
      export_plot(plot.object, plotfilename)
      
      texts_wizard(paste("\nPlot saved in your working directory,",paste(getwd(),", as",sep=""),plotfilename, sep=" "))
    }
    texts_wizard("\nDo you want to try any other PCA plot type?")
    answer <- utils::select.list(choices = truefalse)
  }
  
  
  #### Exporting session log ####
  texts_wizard("\nPCA_analysis_Wizard: Duty Accomplished.")
  dog <- as.character(paste("Session log, ",sinkfilename, " has been saved in your working directory: ", getwd(),sep=""))
  texts_wizard(dog)
  sink()
  sink()
}




