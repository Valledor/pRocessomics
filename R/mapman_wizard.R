#
# Copyright 2019 Luis Valledor / Laura Lamelas
# Last revision 19.04.2020
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
#' @name mapman_wizard
#' @title mapman_wizard() easily analyze and depict mapman categories of your data
#' @description This function is a wizard to interactively help the user to perform the mapman analysis on preprocessed and annotated data
#' @usage mapman_wizard(outputobjectname)
#' @param outputobjectname character, name for mapman result object 
#' @seealso \code{\link{importannotation}}, \code{\link{mapman_group}}, \code{\link{mapman_plot}}
#' @return returns a "MMOG" object and its depiction
#' @author Luis Valledor and Laura Lamelas
#' @export
#'


mapman_wizard<-function(outputobjectname=NULL){
  
  #### Starting console logging ####
  diahora <-gsub(":","",gsub(" ","_", gsub("-","",Sys.time())))
  diahora <-substr(diahora,3,nchar(diahora)-2)
  sinkfilename <- paste(diahora,"_MapmanLog.txt",sep="")
  sink(sinkfilename, split = T)
  
  #### Initial texts ####
  headings_wizards("Welcome to pRocessomics Mapman Wizard")
  texts_wizard("\nThis wizard is aimed to guide you for properly analyse mapman categories within your data. Please note that this module will only work if you have preprocessed your data with preprocess_wizard or preprocess_omic_list functions. Please read 'https://github.com/Valledor/pRocessomics/wiki/heatmaps-and-mapman' before continuing.")
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
  tempdatasetname <- datalists
  
  #### Settings ----
  settings<-selectsettings(originaldatasetname = datalists,datalist = datalist)
  initialrow <- settings[1]
  initialcolumn <- settings[2]
  treatment1col <- settings[3]
  treatment2col <- settings[4]
  treatment <- settings[5]
  
  #### Select omic levels----
  datasetnames<-names(datalist)
  if(length(datasetnames)>1){
    texts_wizard("\nWe have detected the following levels in your data. Please enter the number(s) of the level(s) you want to process and press [enter].\n")
    datasetnames <- c(datasetnames)
    omiclevel <- utils::select.list(datasetnames, preselect = NULL, multiple = T, title = NULL, graphics = FALSE)
  } else {
    omiclevel <- datasetnames
  }
  
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
  
  # Loading annotations
  texts_wizard("\nNow you need to choose the annotations you want to employ to build MapMan groups")
  annotation<-selectannotation()
  if(is.null(annotation)) stop("Fatal error. Annotations are required for continuing")
  class(annotation)<-"data.frame"
  annotation<-as.data.frame(annotation)
  
  #### Mapman_group details.... ####
  texts_wizard("\nDefine the minimum number of variables required for considering a MapMan category (default = 1) and press [enter].")
  threshold <- as.numeric(readline())
  
  texts_wizard("\nDo you want to group your results by treatments? (means of all replicates of each treatment)")
  group<-select.list.edit(choices = truefalse,preselect = NULL,multiple = F,title = NULL,graphics = FALSE)
  if(group == "Yes") group <- TRUE
  if(group =="No") group <- FALSE
  if(group==T){
    texts_wizard("\nDo you want to do statistics over your data and get mean, sd, and pvalues?")
    stats<-select.list.edit(choices = truefalse,preselect = NULL,multiple = F,title = NULL,graphics = FALSE)
    if(stats == "Yes") stats <- TRUE
    if(stats =="No") stats <- FALSE
    
    texts_wizard("\nDo you want scale data for doing stats?") #Elaborar mas la explicacion
    scalevalues<-select.list.edit(choices = truefalse,preselect = NULL,multiple = F,title = NULL,graphics = FALSE)
    if(scalevalues =="No") scalevalues <- FALSE
    #normbyvariables <- FALSE
    if(scalevalues == "Yes") {
      scalevalues <- TRUE}
      # texts_wizard("\nDo you want to scale by variable or by sample?") #esto no lo pillo. basicamente es dejar como esta o normalziar por variable, no?
      # normbyvariables<-select.list.edit(choices = c("Variables","Samples"),preselect = NULL,multiple = F,title = NULL,graphics = FALSE)
      # if(normbyvariables == "Variables") normbyvariables <- TRUE
      # if(normbyvariables =="Samples") normbyvariables <- FALSE}
  } 
  
  #### Call mapman_group ####
  texts_wizard("\nLet's classify your data!")
  
  mapman_grouped_data <- mapman_group(datalist=datalist, annotation = annotation, initialrow=initialrow, initialcolumn=initialcolumn, 
                                      treatment1col=treatment1col,treatment2col=treatment2col,treatment=treatment,
                                      threshold = threshold, group=group, stats = stats, omiclevel = omiclevel, scalevalues = scalevalues)#, normbyvariables = normbyvariables)
  texts_wizard("Press [enter] to continue")
  invisible(readline())
  
  
  #### Naming outputs ####
  if(is.null(outputobjectname)==TRUE) {
    outputobjectname <- paste(datalists,"_mapman",sep="")
    dog<-as.character(paste("\n\nOutput object name was not set by user. Employing default name: ",paste(datalists,"_mapman",sep="")))
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
    outputobjectname <- paste(datalists,"_mapman",sep="")
    dog<-as.character(paste("\n\nOutput object name was wrong. Employing default name: ",paste(tempdatasetname,"_mapman",sep="")))
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
  exportlist <- list(mapman_grouped_data)
  names(exportlist) <- c(outputobjectname)
  list2env(exportlist,envir = .GlobalEnv)
  dog <- as.character(paste(outputobjectname,"has been exported successfuly to your global environment",sep=" "))
  texts_wizard(dog)
  
  texts_wizard("\nDo you want to save MapMan results in .xlsx?")
  answer <- utils::select.list(choices = truefalse)
  if(answer == "Yes"){
    texts_wizard("Saving...")
    diahora <-gsub(":","",gsub(" ","_", gsub("-","",Sys.time())))
    diahora <-substr(diahora,3,nchar(diahora)-2)
    tablefilename<-paste(diahora,"_mapman_results.xlsx",sep="")
    export_table(result_list = mapman_grouped_data, filename = tablefilename)
    texts_wizard(paste("\nResults saved in your working directory,",paste(getwd(),", as",sep=""), tablefilename,".", sep=" "))
  }
  
  
  #### Guided plotting ####
  texts_wizard("\nDo you want to plot your mapman analysis?")
  answer <- utils::select.list(choices = truefalse)
  while (answer == "Yes"){
    texts_wizard("\nWhich type of plot do you want to try?")
    plotchoices <- c("bar","pie","heatmap")
    plottype <- utils::select.list(choices = plotchoices)
    if(plottype == "pie"){
      plot<-mapman_plot(mapman_grouped_data, plottype = "pie",fontsizes = fontsizes)
    }
    if(plottype == "bar"){
      texts_wizard("\nDo you want to logscale Y axis?")
      logscales <- utils::select.list(choices = truefalse)
      if(logscales == "Yes") logscale <- TRUE
      if(logscales == "No") logscale <- FALSE
      texts_wizard("\nDo you want to add error bars?")
      ploterror <- utils::select.list(choices = truefalse)
      if(ploterror == "Yes") ploterror <- TRUE
      if(ploterror == "No") ploterror <- FALSE
      plot<-mapman_plot(mapman_grouped_data, plottype = "bar",logscale = logscale, ploterror = ploterror,fontsizes = fontsizes)
    }
    if(plottype =="heatmap"){
      texts_wizard("\nWhich method de you want to emply for calculating the distance in rows and cols dendrogram in the heatmap?")
      heatmapchoices <- c("correlation","euclidean","maximum","manhattan","camberra","binary","minkowski")
      hm1 <- utils::select.list(choices = heatmapchoices)
      texts_wizard("\nWhich method de you want to emply for sample/treatment clustering?")
      heatmapchoices<-c("ward.D","ward.D2","single","complete","average","mcquiry","median","centroid")
      hm3 <- utils::select.list(choices = heatmapchoices)
      texts_wizard("\nDo you want to display the numbers inside the heatmap cells?")
      answer <- utils::select.list(choices = truefalse)
      if(answer == "Yes") hm4 <- TRUE
      if(answer == "No") hm4 <- FALSE
      texts_wizard("\nWhat about the direction of centering and scaling'")
      choices<-c("row", "column","none")
      hm5 <- utils::select.list(choices = choices)
      hmparameters <- c(hm1,hm1,hm3,hm4,hm5) #hm1 is duplicated in purpose, not a mistake
      plot <- mapman_plot(mapman_grouped_data,plottype = "heatmap", hmparameters = hmparameters,fontsizes = fontsizes)
    }
    print(plot)
    texts_wizard("\nDo you want to save the plot?")
    answer <- utils::select.list(choices = truefalse)
    if(answer == "Yes"){
      texts_wizard("Saving...")
      diahora <-gsub(":","",gsub(" ","_", gsub("-","",Sys.time())))
      diahora <-substr(diahora,3,nchar(diahora)-2)
      plotfilename <- as.character(paste(diahora,plottype,"_mapman_plot",sep=""))
      export_plot(plot_object = plot, filename = plotfilename)
      texts_wizard(paste("\nPlot saved in your working directory,",paste(getwd(),", as",sep=""),plotfilename,".", sep=" "))
      
    }
    texts_wizard("\nDo you want to draw any other plot?")
    answer <- utils::select.list(choices = truefalse)
  }
  #### Exporting session log ####
  texts_wizard("\nMapman_Wizard: Duty Accomplished.")
  dog <- as.character(paste("Session log, ",sinkfilename, " has been saved in your working directory: ", getwd(),sep=""))
  texts_wizard(dog)
  sink()
  sink()
}