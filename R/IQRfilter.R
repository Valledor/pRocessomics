#
# Copyright 2019 Luis Valledor / Laura Lamelas
# Last revision 4.11.2019
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

#' @importFrom stats IQR

IQRfilter <- function(matriz, initialrow, initialcolumn, treatment1col, treatment2col, treatment=1, threshold=0.25,abovebelow="Below"){
  
  if(treatment!=3&treatment!=0){
    if(treatment==1) treatmentcolumn <-treatment1col
    if(treatment==2) treatmentcolumn <-treatment2col
    if(length(unique(matriz[,treatmentcolumn]))==1) stop("IQRfilter Error: Data with more than one level is required")
    datos <- split(as.data.frame(matriz),matriz[,treatmentcolumn],drop=T)
    mediaportratamiento <-lapply(datos,function(x) colMeans(data.matrix(x[initialrow:nrow(x),initialcolumn:ncol(x)]),na.rm=TRUE)) #Sacamos los valores medios de cada tratamiento
    mediaportratamiento <- do.call(rbind,mediaportratamiento) 
    
  }
  if(treatment==3){
    splitvector <- as.vector(paste(matriz[initialrow:nrow(matriz),treatment1col], matriz[initialrow:nrow(matriz),treatment2col],sep="&"))
    datos <- split(as.data.frame(matriz),splitvector,drop=T)
    mediaportratamiento <-lapply(datos,function(x) colMeans(data.matrix(x[initialrow:nrow(x),initialcolumn:ncol(x)]),na.rm=TRUE)) #Sacamos los valores medios de cada tratamiento
    mediaportratamiento <- do.call(rbind,mediaportratamiento) 
  }
  if(treatment==0){ 
    mediaportratamiento <-data.matrix(matriz[initialrow:nrow(matriz),initialcolumn:ncol(matriz)])
  }
  
  iqr <- apply(mediaportratamiento,2,function(x) stats::IQR(x,type = 8)/mean(x)) #
  if(abovebelow=="Above") {
    seleccion<- iqr >= mean(iqr)*(1+threshold) 
  }
  if(abovebelow=="Below"){
    seleccion<- iqr <= mean(iqr)*(1+threshold) 
  }
 
  dog <-paste("\nIQR filtering: Out of",length(seleccion),"variables,",length(seleccion[seleccion==T]),"showed an IQR",abovebelow, threshold*100,"% average IQR.",sep=" ")
  texts_wizard(dog)
  
  matrizfiltrada <- cbind(matriz[,1:initialcolumn-1], matriz[,initialcolumn:ncol(matriz)][,seleccion]) 
  colnames(matrizfiltrada)[1:initialcolumn-1]<-colnames(matriz)[1:initialcolumn-1]
  return(matrizfiltrada)
}
