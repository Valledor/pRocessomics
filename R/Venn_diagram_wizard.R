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

#' Plot Venn diagrams
#' @return The Venn diagram table and plot of an omic dataset
#' @author Luis Valledor and Laura Lamelas
#' @export
#' @importFrom utils select.list
#' @importFrom grDevices pdf

Venn_diagram_wizard<-function(){
  #### Starting console logging ####
  diahora <-gsub(":","",gsub(" ","_", gsub("-","",Sys.time())))
  diahora <-substr(diahora,3,nchar(diahora)-2)
  sinkfilename <- paste(diahora,"_VennWzrdLog.txt",sep="")
  sink(sinkfilename, split = T)
  
  #### Initial texts ####
  headings_wizards("Welcome to pRocessomics Venn Diagram Wizard")
  texts_wizard("\n\nThis wizard is aimed to guide you for ploting Venn Diagrams and generate information about qualitative differences found in your dataset(s). Please read 'https://github.com/Valledor/pRocessomics/wiki/venn' before continuing")
  texts_wizard("\nIf you have preprocessed the data by yourself please change the class of your datalist to 'POL' before proceeding. But do it at your own risk, no data structure and quality check will be performed. We strongly recommend using pRocessomics import, annotation, preprocess and transform functions.")
  texts_wizard("\nPress [enter] to continue...")
  invisible(readline())
 
  #### Select datalist ####
  # We ask the user to input one datalist to be used
  datalists <- Filter(function(x) "POL" %in% class(get(x)), ls(envir = .GlobalEnv))
  if(length(datalists)==0) stop("No valid POL objects were found. Please import your data with importfromexcel().")
  texts_wizard("\nWe have detected the following datasets available. Please enter the number of the dataset do you want to process and press [enter]. If your dataset is not listed please exit pressing [esc] and re-import your data.")
  datalists<- utils::select.list(datalists, preselect = NULL, multiple = F, title = NULL, graphics = FALSE)
  datalist <-get(datalists,envir=.GlobalEnv)
  dog<-as.character(paste("You have selected",datalists,"dataset. This dataset has the following omic levels:\n",sep=" "))
  texts_wizard(dog)
  cat(names(datalist))
  
  #### Settings ----
  settings<-selectsettings(originaldatasetname = datalists,datalist = datalist)
  initialrow <- settings[1]
  initialcolumn <- settings[2]
  treatment1col <- settings[3]
  treatment2col <- settings[4]
  
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
  
  ### Omic levels and so on ----     
  texts_wizard("\nPlease select the omic level(s) you want to include in analysis:")
  datasetnames<-names(datalist)
  omiclevel <- utils::select.list(datasetnames, preselect = NULL, multiple = T, title = NULL, graphics = FALSE)

  
  # Treatment definition and replicates
  if(is.na(treatment2col)){
    treatment<-1
    vectortratamientos <- datalist[[omiclevel[1]]][initialrow:nrow(datalist[[omiclevel[1]]]),treatment1col]
    vectortratamientos <- factor(vectortratamientos, levels = unique(vectortratamientos))
    }
  if(!is.na(treatment2col)){
    texts_wizard("\nPlease select which treatment you want to use for splitting your data.\n")
    treatment<-select.list.edit(c("Treatment 1", "Treatment 2","Combination of Treatments 1 and 2"), preselect = NULL, multiple = F, title = NULL, graphics = FALSE)
    if(treatment=="Treatment 1") {
      treatment<-1
      vectortratamientos <- datalist[[omiclevel[1]]][initialrow:nrow(datalist[[omiclevel[1]]]),treatment1col]
      vectortratamientos <- factor(vectortratamientos, levels = unique(vectortratamientos))
      
    }
    
    if(treatment=="Treatment 2") {
      treatment<-2
      vectortratamientos <- datalist[[omiclevel[1]]][initialrow:nrow(datalist[[omiclevel[1]]]),treatment2col]
      vectortratamientos <- factor(vectortratamientos, levels = unique(vectortratamientos))
    }
    if(treatment=="Combination of Treatments 1 and 2") {
      treatment <-3
      vectortratamientos <- paste(datalist[[omiclevel[1]]][initialrow:nrow(datalist[[omiclevel[1]]]),treatment1col],datalist[[omiclevel[1]]][initialrow:nrow(datalist[[omiclevel[1]]]),treatment2col], sep="&")
      vectortratamientos <- factor(vectortratamientos, levels = unique(vectortratamientos))
    }
  }
 
  replicates <- length(which(vectortratamientos==levels(vectortratamientos)[1]))
  
  ####Threshold ----
  texts_wizard(paste("\nPlease select presence threshold. Select the minimum number of non-zero values per treatment to consider this variable as present. Since you have ", replicates, " replicates please introduce a number between 1 and ", replicates, " and press [enter].",sep=""))
  threshold<-as.numeric(readline())
  
  if(threshold >replicates | threshold < 1){ stop("FATAL ERROR: Please introduce a correct threshold")}
  
  vennmatrix<-Venn_analysis(datalist=datalist,initialrow=initialrow, initialcolumn=initialcolumn, treatment1col=treatment1col,treatment2col=treatment2col, treatment=treatment,omiclevel=omiclevel,threshold=threshold)
  
  a<-Venn_plot(datalist=vennmatrix, Euler.dist = FALSE, alpha=0.5, num.cex=3,font.cex=3)
  
  outputname<-paste(datalists,"_VennGroups", sep="")
  texts_wizard(paste("\n\nA Venn plot summarizing your data has been plot. The different groups employed to draw this plot will be available in ", outputname, ". Press [enter] to continue.", sep=""))
  invisible(readline())
  
  #### Exporting to global environment ------
  dog <- as.character(paste("\n",outputname, " has been saved in global environment.", sep=""))
  exportlist<-list(vennmatrix)
  names(exportlist)<-c(outputname)
  list2env(exportlist,envir = .GlobalEnv)
  texts_wizard(dog)
  # b<- list(a)
  # names(b) <- c("Venn_plot") 
  # list2env(b,envir = .GlobalEnv)
  
  #### Exporting ####
  texts_wizard("\nDo you want to save and export tables and plots?")
  answer <- utils::select.list(choices = truefalse)
  if(answer == "Yes"){
    
    texts_wizard("Saving...")
    diahora <-gsub(":","",gsub(" ","_", gsub("-","",Sys.time())))
    diahora <-substr(diahora,3,nchar(diahora)-2)
    plotfilename <- as.character(paste(diahora,"_venn_plot.pdf",sep=""))
    if(inherits(a,"Venngrid")) export_plot(a,plotfilename)
    else{ #no le quiero meter mas argumentos a export_plot y por eso lo pongo aqui
      grDevices::pdf(plotfilename)
      invisible(a<-Venn_plot(datalist=vennmatrix, Euler.dist = FALSE, alpha=0.5, num.cex=3,font.cex=3))
      dev.off()
    }
    tablefilename<-paste(diahora,"_venn_results.xlsx",sep="")
    export_table(result_list = vennmatrix, filename = tablefilename)
    texts_wizard(paste("\nPlot and table saved in your working directory,",paste(getwd(),", as",sep=""),plotfilename,"and", tablefilename,".", sep=" "))
    
  }
  
 
  #### Exporting session log ####
  texts_wizard("\nVenn_Wizard: Duty Accomplished.")
  dog <- as.character(paste("Session log, ",sinkfilename, " has been saved in your working directory: ", getwd(),sep=""))
  texts_wizard(dog)
  sink()
  sink()
}
