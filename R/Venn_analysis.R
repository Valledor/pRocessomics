#' @name Venn_analysis
#' @title Venn analysis function
#' @description A function to sort variables in order of their presence or absence within the datasets
#' @usage Venn_analysis(datalist,initialrow=1, initialcolumn=2, 
#' treatment1col=1,treatment2col=NA, treatment=1,omiclevel=NULL,threshold=2)
#' @param datalist List with different preprocessed omic levels. pRoDS class object.
#' @param initialrow First row of numeric data within datasets.
#' @param initialcolumn First column of numeric data within datasets.
#' @param treatment1col Column number in which treatment 1 is indicated.
#' @param treatment2col Column number in which treatment 2 is indicated (if applicable).
#' If there is only one treatment indicate the same column as treatment1col.
#' @param treatment numeric. Set 1 to split the dataset according to treatment1, 2 to split the
#' dataset according to treatment 2 and 3 to split the dataset according to the combination of both treatments.
#' @param omiclevel Vector containing the name/s of omics layers to be analyzed
#' @param threshold Numeric indicating the minimum cases in which a variable has to be different to 0 to be taken into account
#' @return A "Vennanalysis" class list object containing the Venn table
#' @author Luis Valledor and Laura Lamelas
#' 



#' @importFrom plyr dlply count
#' @export

Venn_analysis<-function(datalist,initialrow=1, initialcolumn=2, treatment1col=1,treatment2col=NA, treatment=1,omiclevel=NULL,threshold=2){
  #requireNamespace(plyr,quietly = T)
  options(warn=-1)
  if(!inherits(datalist,"POL")) stop("A POL object is required. Please run pre-processing steps")
  datasetnames<-names(datalist)
  
  #Generamos un vector para hacer el split por tratamientos
  if(treatment==1) vectortratamientos <- factor(datalist[[1]][initialrow:nrow(datalist[[1]]),treatment1col],levels = unique(datalist[[1]][initialrow:nrow(datalist[[1]]),treatment1col]))
  if(treatment==2) vectortratamientos <- factor(datalist[[1]][initialrow:nrow(datalist[[1]]),treatment2col],levels = unique(datalist[[1]][initialrow:nrow(datalist[[1]]),treatment2col]))
  if(treatment==3) vectortratamientos <- factor(paste(datalist[[1]][initialrow:nrow(datalist[[1]]),treatment1col],datalist[[1]][initialrow:nrow(datalist[[1]]),treatment2col], sep="&"), 
                                                levels = unique(paste(datalist[[1]][initialrow:nrow(datalist[[1]]),treatment1col],datalist[[1]][initialrow:nrow(datalist[[1]]),treatment2col], sep="&")))
  
  #Extraemos los valores significativos de las tablas
  datalist2 <- lapply(datalist, function(x) x[initialrow:nrow(x),initialcolumn:ncol(x)])
  if(is.null(omiclevel)) omiclevel="all"

  
  # Seleccionamos los niveles omicos datos
  if(omiclevel=="all") {
    variablenames <- lapply(datalist2, function(x) colnames(x))
    variablenames <-unlist(variablenames)
    listaintermedia <-do.call(cbind,datalist2)  #Esto es para evitar cosas raras.
    colnames(listaintermedia) <- variablenames #Esto es para evitar cosas raras.
    listaintermedia<-RemoveEmptyColumns(listaintermedia,1,1) #Esto es para evitar cosas raras.
  }  
  if(omiclevel!="all") {
    sapply(omiclevel, function(x) if(x %in% datasetnames==FALSE) stop("Error: Please select a valid datasetname"))
    
    if(length(omiclevel)>1) {  
      listaintermedia <-do.call(cbind,datalist2[omiclevel])
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
  
  datos <- split(as.data.frame(listaintermedia),vectortratamientos,drop=F)
  
  Treat.data <- lapply(datos,function(x) apply(x,2,function(y) length(which(y != 0))))
  
  Treat.data <- Treat.data[unique(vectortratamientos)]
  presence.matrix <- (do.call("cbind", Treat.data)) 
  
  #hacemos una matriz con los nombres de los tratamientos que estan presentes, t.q. xij>=threshold-->nombredeltreatment; xij<threshold-->NA
  for(i in 1:nrow(presence.matrix)){
    for(j in 1:ncol(presence.matrix)){
      if(presence.matrix[i,j]>=threshold){
        presence.matrix[i,j]<-colnames(presence.matrix)[j]
      } 
      else {
        presence.matrix[i,j]<-NA
      }
    }
  }  
  presence.matrix<-as.data.frame(presence.matrix)

  
  if (length(which(apply(presence.matrix,1,function(x)all(is.na(x)))))>0) {
    cat("\n",length(which(apply(presence.matrix,1,function(x)all(is.na(x))))),"variables have been removed from the analysis, as were not in at least of",threshold,"replicates in any treatment.")
    presence.matrix<-presence.matrix[-which(apply(presence.matrix,1,function(x)all(is.na(x)))),]
  }else{
    cat("\nAll the variales were found to be in at least",threshold,"replicates in one or more treatments.")
  }

 
  grouping <- as.data.frame(apply(presence.matrix,1,function(x) venn_auxiliarpaste(x)))
  colnames(grouping)<-"groups"

  grouping<-cbind(row.names(grouping),grouping) #sip, esto es una cutrez, pero lo necesito para que dlply no me descojone la tabla y luego no haya forma de saber quien es quien
  groups<-colnames(grouping)


  Vennmatrix<-(plyr::dlply(grouping, "groups")) #al split especial este no le sigue un unique porque me da igual el orden


  
  Vennmatrix<-do.call(rbind,Vennmatrix) #de hecho esta mejor ordenado por R que como estaba inicialmente :)
  
  #hasta este punto tengo una matriz de tantas observaciones como variables y cada una con su grupo asociado
  row.names(Vennmatrix)<-Vennmatrix[,1]
  Vennmatrix <- Vennmatrix[,-1]
  
  finaltable<-as.data.frame(plyr::count(Vennmatrix))
  #finaltable nos da los grupos que hay y su frecuencia....
  #ahora lo suyo seria emparejar las variables con sus respectivos grupos...
  finallist<-sapply(finaltable[,1],function(x) row.names(grouping[grouping$groups==x,]))
  
  finallist <- lapply(finallist,function(x) paste(x,collapse="; "))
  Vennfinalobject<-do.call(rbind,finallist)
  Vennfinalobject<-cbind(finaltable,Vennfinalobject)
  colnames(Vennfinalobject)<-c("Groups","Frecuency","Variables")
  #y asi consigo una matriz con el nombre de los grupos de Venn, que luego me podra reconocer el plot, su frecuencia y las variables que lo componen
  if(length(unique(Vennfinalobject$Groups))==1){
    a<-paste(Vennfinalobject$Variables,collapse = ";")
    Vennfinalobject<-as.data.frame(cbind(as.character(Vennfinalobject$Groups[1]),Vennfinalobject$Frecuency[1],a))
    colnames(Vennfinalobject)<-c("Groups","Frecuency","Variables")
  }
  #View(Vennfinalobject)
  vennresult<-list(treatments=vectortratamientos,data=Vennfinalobject)
  class(vennresult) <- "vennanalysis"
  return(vennresult)
}

# Cambiar nombre de argumentos num.cex font.cex para homogeneizarlos con el resto de funciones graficas
