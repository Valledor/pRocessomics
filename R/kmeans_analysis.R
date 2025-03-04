

#' @importFrom dplyr left_join
#' @importFrom plotly plot_ly layout subplot
#' @export
#' 
#' @name kmeans_analysis
#' @title Kmeans analysis
#' @description A function to perform a k-means analysis
#' @usage kmeans_analysis(datalist,annotation=NULL,initialrow=1, initialcolumn=3, 
#' treatment1col=1, treatment2col=2, treatment=1,omiclevel="all",scalation=2,
#' clusters=c(10,20),show.elbow.plot=FALSE)
#' @param datalist List with different preprocessed omic levels. POL class object.
#' @param annotation Character, unquoted name of annotations object, an object generated by importannotation() function of class "annot"
#' @param initialrow First row of numeric data within datasets.
#' @param initialcolumn First column of numeric data within datasets.
#' @param treatment1col Column number in which treatment 1 is indicated.
#' @param treatment2col Column number in which treatment 2 is indicated (if applicable).
#' If there is only one treatment indicate the same column as treatment1col.
#' @param treatment numeric. Set 1 to split the dataset according to treatment1, 2 to split the
#' dataset according to treatment 2 and 3 to split the dataset according to the combination of both treatments.
#' @param omiclevel Character vector indicating the quoted name or names of the omiclevel(s) to be analyzed
#' @param scalation Numeric, indicating the type of scaling method to be applied to the data
#' Available scaling methods for k-means analysis are:
#' scaling but not centering  <- 1
#' scaling and centering <- 2
#' column scaling expressed as decimal <- 3
#' row scaling expressed as decimal <- 4
#' @param clusters Vector (length 2) indicating the maximum and the minimum number of cluster to test
#' @param show.elbow.plot Boolean indicating whether elbow plot is desired to be displayed
#' @return A "kmeananalysis" object, a list with the output of the kmeans analysis
#' @author Luis Valledor and Laura Lamelas



kmeans_analysis<-function(datalist,annotation=NULL,initialrow=1, initialcolumn=3, treatment1col=1, treatment2col=2, treatment=1,omiclevel="all",scalation=2,clusters=c(10,20),show.elbow.plot=FALSE){
  #require(dplyr)
  options(warn=-1)
  if(!inherits(datalist, "POL")) stop("A POL object is required. Please run pre-processing steps")
  datasetnames<-names(datalist)
  
  # split by treatment
  #####ACHTUNG!!!!!!!!! lau FLAG recheck!!!!
  if(treatment==1) vectortratamientos <- factor(datalist[[1]][initialrow:nrow(datalist[[1]]),treatment1col], levels(unique(datalist[[1]][initialrow:nrow(datalist[[1]]),treatment1col])))
  if(treatment==2) vectortratamientos <- factor(datalist[[1]][initialrow:nrow(datalist[[1]]),treatment2col], levels(unique(datalist[[1]][initialrow:nrow(datalist[[1]]),treatment2col])))
  if(treatment==3) vectortratamientos <- factor(paste(datalist[[1]][initialrow:nrow(datalist[[1]]),treatment1col],datalist[[1]][initialrow:nrow(datalist[[1]]),treatment2col], sep="&"), levels = unique(paste(datalist[[1]][initialrow:nrow(datalist[[1]]),treatment1col],datalist[[1]][initialrow:nrow(datalist[[1]]),treatment2col], sep="&")))
  #if(treatment==4) vectortratamientos <- as.factor(paste(datalist[[1]][initialrow:nrow(datalist[[1]]),treatment1col],datalist[[1]][initialrow:nrow(datalist[[1]]),treatment2col], sep="&"))
  
 
  #### Treatment table ####
  treatments <- buildtreatments(datalist,c(initialrow, initialcolumn, treatment1col, treatment2col, treatment))
  
  # Numeric values extraction
  datalist2 <- lapply(datalist, function(x) x[initialrow:nrow(x),initialcolumn:ncol(x)])
  if(is.null(omiclevel)==TRUE) omiclevel="all"
  
  # Omic levels selection
  if(omiclevel=="all") {
    variablenames <- lapply(datalist2, function(x) colnames(x))
    variablenames <-unlist(variablenames)
    listaintermedia <-do.call(cbind,datalist2)  #### Esto es inconsistente. Cambia nombres de columnas, a veces, nombretabla.nombrevar, otras veces solo nombrevar
    colnames(listaintermedia) <- variablenames #### para arreglar lode arriba entre tanto.
    listaintermedia<-RemoveEmptyColumns(listaintermedia,1,1)#Esto es para evitar cosas raras.
  }  
  if(omiclevel!="all") {
    sapply(omiclevel, function(x) if(x %in% datasetnames==FALSE) stop("Error: Please select a valid datasetname"))
    
    if(length(omiclevel)>1){  listaintermedia <-do.call(cbind,datalist2[omiclevel])
    variablenames <- lapply(datalist2[omiclevel], function(x) colnames(x))
    variablenames <-unlist(variablenames)
    colnames(listaintermedia) <- variablenames
    listaintermedia<-RemoveEmptyColumns(listaintermedia,1,1)#Esto es para evitar cosas raras.
    }
    if(length(omiclevel)==1){
      listaintermedia <-datalist2[[omiclevel]]
      variablenames <-colnames(listaintermedia)
      listaintermedia<-RemoveEmptyColumns(listaintermedia,1,1)#Esto es para evitar cosas raras.
    }
  }
  
  if(!is.null(annotation)){ #Annotation table
    #library(dplyr)
    annotation <- annotation[,1:4]
    colnames(annotation) <- c("IDENTIFIER","DESCRIPTION","MAPMAN.BINCODE","MAPMAN.BINNAME")
    columnnames <- c("IDENTIFIER","DESCRIPTION","MAPMAN.BINCODE","MAPMAN.BINNAME")
    annotmatrix <-cbind(variablenames,"") # required to keep table structure when using dplyr
    colnames(annotmatrix) <- c("IDENTIFIER","")
    annotmatrix <- dplyr::left_join(as.data.frame(annotmatrix),annotation,by="IDENTIFIER") #change to merge!
    annotmatrix <- annotmatrix[,-2]
    annotmatrix <- sapply(annotmatrix,as.vector)
    annotmatrix[annotmatrix[,2]=="",2]<-"" #if no annot <- empty
    annotmatrix[annotmatrix[,3]=="",3]<- 35.2 #if NA ... unknown anyway
    annotmatrix[annotmatrix[,4]=="",4]<-"not assigned" # unknown anyway
    annotmatrix[,3] <- as.numeric(do.call(rbind,lapply(strsplit(as.vector(annotmatrix[,3]),"\\."), function(x) x[[1]]))) #first bin
    annotmatrix[,4] <- as.character(do.call(rbind,lapply(strsplit(as.vector(annotmatrix[,4]),"\\."), function(x) x[[1]]))) #description
    colnames(annotmatrix) <- columnnames
    annotmatrix[,2][is.na(annotmatrix[,2])] <-""
    annotmatrix[,3][is.na(annotmatrix[,3])] <-"35"
    annotmatrix[,4][is.na(annotmatrix[,4])] <-"not assigned"
    annotation<-annotmatrix
  }
  datos <- split(as.data.frame(listaintermedia),vectortratamientos,drop=T)
  datos<-datos[unique(vectortratamientos)]
  #a veces en R al editar o generar una matriz los valores numericos se convierten en factores, para volver a valores numericos
  #los defines primero como caracteres y luego como numero, si pones as.numeric  a pelo MAL MuY MAL, porque te devuelve el ordinal 
  #de cada factor, es decir te los ordena de menor a mayor y el numero que te devuelve es la posicion en la lista
  #no da error, y no te das cuenta...
  datos<-lapply(datos, function(datos) t(apply(datos,1,function(datos) as.numeric(as.character(datos)))))
  #media de las replicas por cada variable
  treat.means<-lapply(datos, colMeans)
  #juntarlo en una tabla, gracias al unique no se van a descolocar nombres
  means.matrix<-do.call("rbind",treat.means)
  #y les volvemos a poner los nombres a las variables 
  colnames(means.matrix)<-colnames(listaintermedia)
  #hasta aqui la matriz para el k-means clustering
  #dim(means.matrix)=n? treatments x n? variables
  #Different choices to scale/center/percentage... whatever
  #scale
  if(scalation==1){meansfmatrix<-scale(means.matrix,center=F)}
  #scale and center
  else if (scalation==2){meansfmatrix<-scale(means.matrix)}
  #percentaje (tanto por 1) por columna, a mi me gusta este!
  else if (scalation==3){meansfmatrix<-sweep(means.matrix,2,colSums(means.matrix),"/" )  }
  #percentaje (tanto por 1) por fila, otra opcion
  else if (scalation==4){meansfmatrix<-sweep(means.matrix,1,rowSums(means.matrix),"/" )  }
  
  ## Ahora que los datos estan bien nos ponemos con el kmeans
  groups<-c()
  groups.ensayo <-c()
  elbow.data<-c()
  lista.kmeans <-list()
  kmeans_matrix <- t(meansfmatrix)
  #if(clusters=="auto"){} #sin implementar aqui. La idea es rango de kmeans entre 5 y 50, y seleccionar rango optimo
  # hacemos el analisis kmeans
  j=1
  for(i in (clusters[1]:clusters[2])){
    kmeans.result <-kmeans_helper(kmeans_matrix, i)
    groups.ensayo <- cbind(groups.ensayo,kmeans.result[[1]])
    colnames(groups.ensayo)[j] <- paste(i,"Clusters",sep=" ")
    elbow.data <- c(elbow.data,kmeans.result[[2]])
    names(elbow.data)[j] <-paste(i,"Clusters",sep=" ")
    lista.kmeans[[j]] <-as.data.frame(kmeans.result[[3]])
    names(lista.kmeans)[[j]] <- paste(i,"Clusters",sep=" ")
    j=j+1
  }
  if(show.elbow.plot==TRUE){
    
    toplotly<-as.data.frame(cbind("clusters"=clusters[1]:clusters[2],"Elbow"=elbow.data))
    #View(toplotly)
    p <- elbowplot(toplotly)
   
    
  } 
  if(!is.null(annotation)){
    lista.kmeans <- lapply(lista.kmeans, function(x) {
      merge(x,annotmatrix, by.x="ID", by.y="IDENTIFIER")
    })
    resultado <-list(meansfmatrix,annotmatrix,lista.kmeans,groups.ensayo,elbow.data, treatments)
    names(resultado) <- c("kmeans_matrix","annotations","kmeans_list","kmeans_group","elbow.data", "treatments")}
  else{resultado <-list(meansfmatrix,lista.kmeans,groups.ensayo,elbow.data, treatments)
  names(resultado) <- c("kmeans_matrix","kmeans_list","kmeans_group","elbow.data", "treatments")}
  class(resultado) <- "kmeansanalysis"
  return(resultado)
}
