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

#' @importFrom stats sd

VarCoefSelect <- function(matriz, initialrow, initialcolumn, treatment1col, treatment2col, treatment=1, threshold=0.75,datasetname=NULL,abovebelow="Below"){
  options(warn=-1) 
  
  if(treatment==1) vectortratamientos <- (matriz[initialrow:nrow(matriz),treatment1col])
  if(treatment==2) vectortratamientos <- (matriz[initialrow:nrow(matriz),treatment2col])
  if(treatment==3) vectortratamientos <- (paste(matriz[initialrow:nrow(matriz),treatment1col],matriz[initialrow:nrow(matriz),treatment2col], sep="&"))
  vectortratamientos<-factor(vectortratamientos, levels = unique(vectortratamientos))
  
  datos <- split(as.data.frame(matriz),vectortratamientos,drop=T) 
  coefvariacion <-lapply(datos,function(x) apply(x[initialrow:nrow(x),initialcolumn:ncol(x)],2,function(y) stats::sd(y)/mean(y))) 
  coefvariacion <- do.call("rbind",coefvariacion) 
  coefvariacion <- apply(coefvariacion,2,min,na.rm=T) 
  
  if(abovebelow=="Above") {
    seleccion <- coefvariacion >= threshold
  }
  if(abovebelow=="Below"){
    seleccion <- coefvariacion <= threshold
  }
  
  dog <-paste("\nVarCoef filtering: Out of",length(seleccion),"variables,",length(seleccion[seleccion==T]),"showed a CV", abovebelow, threshold*100,"% in at least one treatment/category.",sep=" ")
  texts_wizard(dog)
  
  matrizfiltrada <- cbind(matriz[,1:initialcolumn-1], matriz[,initialcolumn:ncol(matriz)][,seleccion])
  colnames(matrizfiltrada)[1:initialcolumn-1]<-colnames(matriz)[1:initialcolumn-1]
  return(matrizfiltrada)
}
