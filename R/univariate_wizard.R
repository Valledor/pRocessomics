#
# Copyright 2019 Luis Valledor / Laura Lamelas
# Last revision 08.01.2019
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
#' @name univariate_wizard
#' @title Univariate analysis wizard
#' @description A function to perform univariate analysis within omic datasets
#' @usage univariate_wizard()
#' @details This wizard is aimed to guide you for properly analyze your datasets. Please note that this module will only work if you have preprocessed your data with preprocess_wizard or preprocess_omic_list functions.
#' During the wizard you will be asked to introduce your data, you will be allowed to check parametricity, transform your data and perform parametric and non-parametric tests to evaluate the significance of variables.
#' @seealso univariate, preprocess_omic_list, preprocess_wizard, transformadata, transformation_wizard
#' @author Luis Valledor and Laura Lamelas
#' @export

univariate_wizard<-function(){
  #### Starting console logging ####
  diahora <-gsub(":","",gsub(" ","_", gsub("-","",Sys.time())))
  diahora <-substr(diahora,3,nchar(diahora)-2)
  sinkfilename <- paste(diahora,"_UnivarWzrdLog.txt",sep="")
  sink(sinkfilename, split = T)
  
  #### Initial text ####
  headings_wizards("Welcome to pRocessomics Univariate Analysis Wizard")
  texts_wizard("\n\nThis wizard is aimed to guide you for properly analyze your datasets. Please note that this module will only work if you have preprocessed your data with preprocess_wizard() or preprocess_omic_list functions().Please read 'https://github.com/Valledor/pRocessomics/wiki/Univariate' before continuing.")
  texts_wizard("\nIf you have preprocessed the data by yourself please change the class of your datalist to 'POL' before proceeding. But do it at your own risk, no data structure and quality check will be performed. We strongly recommend using pRocessomics import, annotation, preprocess and transform functions.")
  texts_wizard("\nPress [enter] to continue...")
  invisible(readline())
  
  #### Select datalist ####
  # We ask the user to input one datalist to be used
  datalists <- Filter(function(x) "POL" %in% class(get(x)), ls(envir = .GlobalEnv))
  if(length(datalists)==0) stop("No valid POL objects were found. Please preprocess your data with preprocess_omic_list() function or preprocess_wizard().")
  texts_wizard("\nWe have detected the following datasets available. Please enter the number of the dataset do you want to transform and press [enter].")
  texts_wizard("\nIf your dataset is not listed please exit pressing [esc] and re-import your data.")
  datalists<- select.list.edit(datalists, preselect = NULL, multiple = F, title = NULL, graphics = FALSE)
  datalist <-get(datalists,envir=.GlobalEnv)
  dog<-as.character(paste("You have selected",datalists,"dataset. This dataset has the following omic levels:\n",sep=" "))
  texts_wizard(dog)
  datalist.to.univar<-datalist
  datalist.to.univar.name<-datalists
  cat(names(datalist))
  
  #### Settings ----
  settings<-selectsettings(originaldatasetname = datalists,datalist = datalist)
  initialrow <- settings[1]
  initialcolumn <- settings[2]
  treatment1col <- settings[3]
  treatment2col <- settings[4]
  #treatment <- settings[5]
  
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
  # Treatment definition 
  if(is.na(treatment2col)){treatment=1}
  if(!is.na(treatment2col)){
    texts_wizard("\nPlease select which treatment you want to use for splitting your data.\n")
    treatment<-select.list.edit(c("Treatment 1", "Treatment 2","Combination of Treatments 1 and 2"), preselect = NULL, multiple = F, title = NULL, graphics = FALSE)
    if(treatment=="Treatment 1") treatment<-1
    if(treatment=="Treatment 2") treatment<-2
    if(treatment=="Combination of Treatments 1 and 2") treatment <-3
  }
  
  
  
  #### Intro 2 --------
  texts_wizard("\n\nDo you want to test is your data is parametric (normal-distributed and homocedastic) before starting? Knowing your data before performing univariate analyses is relevant for deciding the best statistical test to apply.")
  #### Are you parametric? -------
  yesno <-c("Yes","No")
  testparametric<- select.list.edit(yesno, preselect = NULL, multiple = F, title = NULL, graphics = FALSE)
  if(testparametric=="Yes"){
    
    parametricity(datalist=datalist, initialcolumn=initialcolumn, initialrow=initialrow, treatment1col=treatment1col,
                  treatment2col=treatment2col, treatment=treatment,name=datalists)
    #texts_wizard("\nPress [enter] to continue...")
    #1invisible(readline())  
    
    # Do you want to try any transform your data?----
    texts_wizard("\nIf you are interested in using parametric tests but your data does not fulfill normality and/or homocedascity assumptions, you can try again with data transformation.")
    texts_wizard("\nData transformation is a common step in omics workflow, in pRocessomics there are several data transformation methods implemented.")
    texts_wizard("\nDo you want to try any data transformation?")
    datatransform<- select.list.edit(yesno, preselect = NULL, multiple = F, title = NULL, graphics = FALSE)
    
    userschoices <- c()
    while(datatransform=="Yes"){
      dog<-as.character(paste("\nYour omic list contains ",length(datalist)," datasets: ", paste(names(datalist),collapse=" "), sep=""))
      texts_wizard(dog)
      texts_wizard(paste("\nPlease select one transformation method for each dataset. It is expected that you type ",length(names(datalist)), " numbers separated by spaces i.e. \"",paste(c(1:length(names(datalist))),collapse=" "),"\". Follow the order of the datasets shown above. If you don't want to transform a particular dataset use 9 (none).",sep=""))
      tchoices<-c("z","z-centered","Log10","Log10plus1","box-cox-lambda","sqrt","cubic","arcsine","none")
      transf.choose<- select.list.edit(tchoices, preselect = NULL, multiple = T, title = NULL, graphics = FALSE)
      if(length(transf.choose) != length(datalist)){
        texts_wizard(paste("\nERROR: You have not typed",length(datalist),"elements for transformation. Going back to previous step. Press [enter] to continue.",sep= " "))
        invisible(readline())
        next()
      }
      
      #transformation<-transf.choose
      userschoices <- cbind(userschoices,transf.choose)
      
      transformed.datalist<-transformdata(datalist=datalist, initialrow=initialrow, initialcolumn=initialcolumn, 
                                          treatment1col=treatment1col,treatment2col=treatment2col,treatment=treatment,
                                          transf=transf.choose)
      
      texts_wizard("\nNow it's time to re-check parametricity.Press [enter] to continue")
      invisible(readline())
      
      name.t<-paste("transformed ",datalists,sep="")
      parametricity(datalist=transformed.datalist, initialcolumn=initialcolumn, initialrow=initialrow, treatment1col=treatment1col,
                    treatment2col=treatment2col, treatment=treatment, name=name.t)
      
      texts_wizard("\nDo you want to try another transformation?")
      datatransform<- select.list.edit(yesno, preselect = NULL, multiple = F, title = NULL, graphics = FALSE)
      
      if(datatransform=="No"){
        texts_wizard("\nPlease, select one of the transformations you have already tested")
        dog<-as.character(paste("\nDataset order: ", paste(names(datalist),collapse=" "), sep=""))
        texts_wizard(dog)
        tchoices <- do.call(paste, split(userschoices, seq(nrow(userschoices))))
        transf.choose<- select.list.edit(tchoices, preselect = NULL, multiple = T, title = NULL, graphics = FALSE)
        transf.choose <- do.call(as.vector, strsplit(transf.choose, " "))
        
        suppressMessages(output.object<-transformdata(datalist=datalist, initialrow=initialrow, initialcolumn=initialcolumn, 
                                                      treatment1col=treatment1col,treatment2col=treatment2col,treatment=treatment,
                                                      transf=transf.choose))
        
        output.object.name<-paste(datalists,"_transformed", sep="")
        dog<-as.character(paste("\nPlease note, ",output.object.name," has been saved in your global environment."), sep="")
        texts_wizard(dog)
        exportlist<-list(output.object,transf.choose)
        applied.transf<-"transformation"
        names(exportlist)<-c(output.object.name,applied.transf)
        list2env(exportlist,envir = .GlobalEnv)
        datalist.to.univar<-output.object
        datalist.to.univar.name<-output.object.name
      }
    }
    
  }
  #Univariate ------
  
  if(treatment==1) splitvector <- as.vector(datalist.to.univar[[1]][,treatment1col])
  if(treatment==2) splitvector <- as.vector(datalist.to.univar[[1]][,treatment2col])
  if(treatment==3) splitvector <- as.vector(paste(datalist.to.univar[[1]][,treatment1col], datalist[[1]][,treatment2col],sep="&"))
  
  stattest <- c()
  if(length(unique(splitvector))==1) stop("Two or more treatments are required")
  if(length(unique(splitvector))==2) stattest<-"t_or_u"
  if(length(unique(splitvector))>=3) stattest<-"anova_or_kruskal"
  
  
  if(stattest=="t_or_u"){
    texts_wizard("\nYou are comparing two treatments. Two univariate tests are available: t-test (parametric) and U-test (non-parametric).")}
  if(stattest=="anova_or_kruskal"){
    texts_wizard("\nYou are comparing more than two treatments. Two univariate tests are available: ANOVA (parametric) and Kruskal-Wallis (non-parametric).")}
  
  texts_wizard("\nPlease select the test you want to apply")
  choices_uni<-c("Parametric","Non-parametric")
  params<-select.list.edit(choices = choices_uni,preselect = NULL,multiple = F,title = NULL,graphics = FALSE)
  if(params=="Parametric"){parametric=TRUE} else {parametric=FALSE}
  posthoc=FALSE
  if(stattest=="anova_or_kruskal"){
    texts_wizard("\nDo you want to perform a posthoc analysis (Tukey's HSD or Duncan) to compare treatments by pairs?")
    posth<-select.list.edit(yesno,preselect = NULL,multiple = F,title = NULL,graphics = FALSE)
    if(posth=="Yes"){posthoc=TRUE} else {posthoc=FALSE}}
  texts_wizard("\nDo you want to calculate FDR (BH)? ")
  fd<-select.list.edit(yesno,preselect = NULL,multiple = F,title = NULL,graphics = F)
  if(fd=="Yes"){FDR=TRUE} else {FDR=FALSE}
  
  anotlists <- Filter(function(x) "pRoAnnot" %in% class(get(x)), ls(envir = .GlobalEnv))
  if(length(anotlists)==0) {
    annotatefile<-NULL
    annotating <- "No"}
  if(length(anotlists)!=0) {
    texts_wizard("\nDo you want to add the variable \"description\" from your annotation file to the resulting tables?")
    annotating<-select.list.edit(yesno,preselect = NULL,multiple = F,title = NULL,graphics = FALSE)
  }
  if(annotating =="Yes"){
    if(length(anotlists)==1){
      annotatefile <-get(anotlists,envir=.GlobalEnv)
    } else {
      texts_wizard("There are more than one annotations in global environment. Please select the one you want to use.")
      anntt<-select.list.edit(choices = anotlists,preselect = NULL,multiple = F,title = NULL,graphics = FALSE)
      annotatefile <-get(anntt,envir=.GlobalEnv)  
    }
  }
  if(annotating =="No"){annotatefile<-NULL}
  
  texts_wizard("\nStarting Univariate Analysis\n")
  univarfin<-univariate(datalist = datalist.to.univar,initialrow = initialrow, initialcolumn = initialcolumn, treatment1col = treatment1col, treatment2col = treatment2col,treatment = treatment,
                        parametric = parametric, posthoc = posthoc, FDR = FDR, round = 4,  annotatefile = annotatefile)
  
  #### Exporting things to global environment ------
  outputname<-paste(datalist.to.univar.name,"_univariate",sep="")
  dog <- as.character(paste("\n\n",outputname, " has been saved in global environment.", sep=""))
  exportlist<-list(univarfin)
  names(exportlist)<-outputname
  list2env(exportlist, envir = .GlobalEnv)
  texts_wizard(dog)
  ######
  #plotling
  texts_wizard("\nDo you want to plot the pvalues of your univariate's results?")
  swp <- utils::select.list(choices = truefalse)
  if(swp == "Yes"){
    texts_wizard("\nPlease select between distribution or histogram")
    plottype <- utils::select.list(choices = c("distribution","histogram"))
    annotatefilep <- NULL
    if(plottype == "distribution" ){
      if(!is.null(annotatefile)){
        texts_wizard("\nDo you want to show variables' annotations within the plot?")
        antp <- utils::select.list(choices = truefalse)
        if(antp == "Yes") annotatefilep = annotatefile
        else annotatefilep <- NULL
      }
    }
    newplot <- univariate_plot(univarobject=univarfin, plottype=plottype, annotatefile= annotatefilep, savetofile=F, filename=NULL)
    
  }
  #### Exporting ####
  texts_wizard("\nDo you want to save univariate analysis tables and associated plots in XLS/PDF?")
  answer <- utils::select.list(choices = truefalse)
  if(answer == "Yes"){
    
    texts_wizard("Saving...")
    diahora <-gsub(":","",gsub(" ","_", gsub("-","",Sys.time())))
    diahora <-substr(diahora,3,nchar(diahora)-2)
    plotfilename <- as.character(paste(diahora,"_univar_plot",sep=""))
    invisible(univariate_plot(univarobject = univarfin ,plottype=plottype,annotatefile=annotatefilep,savetofile=T,filename=plotfilename))
    ## añadir en exportplot el modo para los univar, y que añada un "_nombredetabla.pdf" al filename que le mandamos.
    #export_plot(plot_object = plot, filename = plotfilename)
    
    #res.parametricity<-list(ks.list,levenes) #_aqui hay que export univarfin para metricity ya tuvo su momento xD
    tablefilename<-paste(diahora,"_univariate_results.xlsx",sep="")
    export_table(result_list = univarfin, filename = tablefilename)
    texts_wizard(paste("\nPlot and table saved in your working directory,",paste(getwd(),", as",sep=""),plotfilename,"and", tablefilename,".", sep=" "))
    
  }
  
  #### Exporting session log ####
  texts_wizard("\nUnivariate_Wizard: Duty Accomplished.")
  dog <- as.character(paste("Session log, ",sinkfilename, " has been saved in your working directory: ", getwd(),sep=""))
  texts_wizard(dog)
  sink()
  
}

