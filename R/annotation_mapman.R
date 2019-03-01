annotation_mapman <- function(datos,annotationmatrix,removeduplicates=T){ #ddd
  library(dplyr)
  if(class(datos)%in%c("list","POL")){
    variablenames <- lapply(datos, function(x) colnames(x))
    variablenames <-unlist(variablenames)
  } else {
    variablenames <-colnames(datos)
  }
  columnnames <- c("IDENTIFIER","DESCRIPTION","MAPMAN.BINCODE","MAPMAN.BINNAME")
  annotmatrix <-cbind(variablenames,"") # esto es una caca, pero sino el dplyr no me funciona
  colnames(annotmatrix) <- c("IDENTIFIER","")
  
  annotmatrix <- left_join(as.data.frame(annotmatrix),annotationmatrix,by="IDENTIFIER")
  annotmatrix <- annotmatrix[,-2]
  annotmatrix <- sapply(annotmatrix,as.vector)
  
  annotmatrix[is.na(annotmatrix)] <-""
  annotmatrix[annotmatrix[,2]=="",2]<-annotmatrix[annotmatrix[,2]=="",1] #Dejamos en blanco todo aquello que no tenga anotaciones
  annotmatrix[annotmatrix[,3]=="",3]<-35.2 #Definimos como desconocidos los valores perdidos
  annotmatrix[annotmatrix[,4]=="",4]<-"not assigned" # lo mismo.
  
  annotmatrix[,3] <- as.numeric(do.call(rbind,lapply(strsplit(as.character(annotmatrix[,3]),"\\."), function(x) x[[1]]))) #me quedo con el primer bin
  annotmatrix[,4] <- as.character(do.call(rbind,lapply(strsplit(as.character(annotmatrix[,4]),"\\."), function(x) x[[1]]))) #me quedo con el primer descriptor
  colnames(annotmatrix) <- columnnames
  
  annotmatrix[,3][annotmatrix[,3]==""] <-"35"
  annotmatrix[,4][annotmatrix[,4]==""] <-"not assigned"
  
  
  # This removes duplicated paired ID and mapman group pairs
  if(removeduplicates==T) annotmatrix<-annotmatrix[!duplicated(annotmatrix[,c(1,3)]),]
  
  annotation<-as.matrix(annotmatrix)
  return(annotation)
}
