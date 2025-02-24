#' @importFrom utils select.list
#' @export
#' @name kmeans_wizard
#' @title kmeans wizard
#' @description A function to ease the kmneans analysis
#' @usage kmeans_wizard()
#' @author Luis Valledor and Laura Lamelas
#' @return kmeans result list and kmeans analysis plot

kmeans_wizard<-function(){
  #### Initial text ####
  headings_wizards("Welcome to pRocessomics K means analysis Wizard")
  texts_wizard("\n\nThis wizard is aimed to guide you for perform a k-means analysis")
  texts_wizard("\nPlease note that this module will only work if you have preprocessed your data with preprocess_wizard or preprocess_omic_list functions.\n")
  
  texts_wizard("\nIf you have preprocessed the data by yourself please change the class of your datalist to 'POL' before proceeding. But do it at your own risk, since no data structure and quality check will be performed. We strongly recommend using pRocessomics import, annotation, preprocess, filter and transform functions.\n\n")
  texts_wizard("\nPress [enter] to continue...")
  invisible(readline())
  
  #### Select datalist ####
  # We ask the user to input one datalist to be used
  datalists <- Filter(function(x) "POL" %in% class(get(x)), ls(envir = .GlobalEnv))
  if(length(datalists)==0) stop("No valid POL objects were found. Please preprocess your data with preprocess_omic_list function or preprocess_wizard.")
  texts_wizard("\n We have detected the following datasets available. Please enter the number of the dataset do you want to transform and press [enter]. If your dataset is not listed please exit pressing [esc] and re-import your data.")
  datalists<- select.list.edit(datalists, preselect = NULL, multiple = F, title = NULL, graphics = FALSE)
  writeLines(datalists, con = stdout(), sep = "\n", useBytes = FALSE) 
  
  datalist <-get(datalists,envir=.GlobalEnv)
  dog <-as.character(cat(paste("You have selected",datalists,"dataset. This dataset has the following omic levels:\n")))
  texts_wizard(dog)
  cat(names(datalist))
  #### Settings ----
  settings<-selectsettings(originaldatasetname = datalists,datalist = datalist)
  # if(exists(paste(datalists,"_dsettings",sep=""), envir = .GlobalEnv)==TRUE){
  #   settings <- get(paste(datalists,"_dsettings",sep=""),envir = .GlobalEnv)
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
  # Treatment definition 
  if(is.null(treatment2col)){treatment=1}
  if(!is.null(treatment2col)|!is.na(treatment2col)){
    texts_wizard("\nPlease select which treatment you want to use for splitting your data.\n")
    treatment<-select.list.edit(c("Treatment 1", "Treatment 2","Combination of Treatments 1 and 2"), preselect = NULL, multiple = F, title = NULL, graphics = FALSE)
    writeLines(treatment, con = stdout(), sep = "\n", useBytes = FALSE) 
    
    if(treatment=="Treatment 1") treatment<-1
    if(treatment=="Treatment 2") treatment<-2
    if(treatment=="Combination of Treatments 1 and 2") treatment <-3
  }
  
  ### kmeans ----     
  texts_wizard("\nK-means clustering is a simple unsupervised learning algorithm that is used to cluster variables which follow the same pattern all across the treatments.")
  texts_wizard("\nIt follows a simple procedure of classifying a dataset into a number of clusters defined by the letter k (K-means), which will be fixed later.")
  texts_wizard("\nAll the variables will be then sorted according to the nearest behaviour pattern.")
  texts_wizard("\nTo ease the clustering, variables should be scaled first, since in this kind of nanalysis is focused in general patterns and not in absolute values.")
  texts_wizard("\nAvailable scaling methods for k means analysis are:")
  cat("\nscaling but not centering \nscaling and centering \ncolumn scaling expressed as decimal \nrow scaling expressed as decimal")
  dog<-(paste("\nYour omic list contains ",length(datalist)," datasets:", sep=""))
  texts_wizard(as.character(dog))
  cat(names(datalist))
  kchoices<-c("scaling and NOT centering","scaling and centering","column scaling as dec","row scaling as dec")
  k.choose<- select.list.edit(kchoices, preselect = NULL, multiple = T, title = NULL, graphics = FALSE)
  writeLines(k.choose, con = stdout(), sep = "\n", useBytes = FALSE) 
  
  if(k.choose=="scaling and NOT centering") scalation<-1
  if(k.choose=="scaling and centering") scalation<-2
  if(k.choose=="column scaling as dec") scalation<-3
  if(k.choose=="row scaling as dec") scalation<-4
  
  
  
  ### Omic levels and so on ----     
  texts_wizard("\nPlease select which omic levels de you want to analyze:\n")
  datasetnames<-c(names(datalist),"all")
  omiclevel <- utils::select.list(datasetnames, preselect = NULL, multiple = T, title = NULL, graphics = FALSE)
  writeLines(omiclevel, con = stdout(), sep = "\n", useBytes = FALSE) 
  
  ### elbow.plot from 1 to 100 ----
  texts_wizard("\nTo set the optimal number of clusters, we are going to use the elbow method. ")
  texts_wizard("\nThis method consist of selecting the number of clusters (x axis) where the cumulative whitin-clusters sum of squares (y axis) value drops. ")
  texts_wizard("\nWe strongly recommend to select a range of clusters centered in the theoric optimun eg(10 - 15 if the optimum is in 12 clusters) ")
  
  kmeansobject<-kmeans_analysis(datalist=datalist,annotation=NULL,initialrow=initialrow, initialcolumn=initialcolumn, treatment1col=treatment1col, treatment2col=treatment2col, 
                  treatment=treatment,omiclevel=omiclevel,scalation=scalation,clusters=c(1,100),show.elbow.plot=TRUE)
  
  #return(kmeansobject)
  ### number of cluster selection -----
  texts_wizard("\nPlease select the maximum cluster number value\n")
  clusmax<-readline()
  writeLines(clusmax, con = stdout(), sep = "\n", useBytes = FALSE) 
  
  clusmax<-as.numeric(clusmax)
  texts_wizard("\nPlease select the minimum cluster number value\n")
  clusmin<-readline()
  writeLines(clusmin, con = stdout(), sep = "\n", useBytes = FALSE) 
  
  clusmin<-as.numeric(clusmin)
  
  cluster<-c(clusmin,clusmax)
  if(cluster[1]==cluster[2]){cluster<-clusmin}
  ### Call the main function ----
  kmeansobject<-kmeans_analysis(datalist=datalist,annotation=NULL,initialrow=initialrow, initialcolumn=initialcolumn, treatment1col=treatment1col, treatment2col=treatment2col, 
                                   treatment=treatment,omiclevel=omiclevel,scalation=scalation,clusters=cluster,show.elbow.plot=FALSE)
  texts_wizard("\nkmeans analysis....Done!\n")

 ###### FLAG FALG FLAG voy por aqui hace falta meter si no detras de la preguna etc etc, y seguir capturando las elecciones del user 
  texts_wizard("\nDo you want to plot your kmeans results?")
  
  #### Plot ----
  
  a <- cluster[1]:cluster[2]
  
  
  for (i in 1:length(a)){
    kmeans_plot(kmeansobject,clusternumber=(a[i]))
  }
  
  
  #### Export and end ---
  output.object.name<-paste(datalists,"_kmeans", sep="")
  dos<-as.character(paste("\nPlease note, ",output.object.name," has been saved in your global environment."), sep="")
  texts_wizard(dos)
  exportlist<-list(kmeansobject)
  names(exportlist)<-c(output.object.name)
  list2env(exportlist,envir = .GlobalEnv)
  
}

