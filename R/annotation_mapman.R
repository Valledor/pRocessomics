#
# Copyright 2019 Luis Valledor / Laura Lamelas
# Last revision 31.05.2024
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

#' @importFrom dplyr left_join

###### try to remove dplyr left_join 
##### re-write to allow for annotate analysis objects
##### cuidado con el if %in%

annotation_mapman <- function(datos,annotationmatrix,removeduplicates=T){ 
  #library(dplyr)
  if(class(datos)%in%c("list","POL")){. 
    variablenames <- lapply(datos, function(x) colnames(x))
    variablenames <-unlist(variablenames)
  } else {
    variablenames <-colnames(datos)
  }
  columnnames <- c("IDENTIFIER","DESCRIPTION","MAPMAN.BINCODE","MAPMAN.BINNAME")
  annotmatrix <-cbind(variablenames,"") # esto es una caca, pero sino el dplyr no me funciona
  colnames(annotmatrix) <- c("IDENTIFIER","")
  
  annotmatrix <- dplyr::left_join(as.data.frame(annotmatrix),annotationmatrix,by="IDENTIFIER")
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
