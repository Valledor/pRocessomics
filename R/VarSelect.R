#
# Copyright 2019 Luis Valledor / Laura Lamelas
# Last revision 04.06.2019
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


VarSelect <- function(matriz, initialrow, initialcolumn, treatment1col, treatment2col, treatment=1, threshold=0.25,datasetname=NULL){
  cat(paste("\nProcessing ",datasetname," dataset",sep=""))
  if(treatment==1) vectortratamientos <- (matriz[initialrow:nrow(matriz),treatment1col])
  if(treatment==2) vectortratamientos <- (matriz[initialrow:nrow(matriz),treatment2col])
  if(treatment==3) vectortratamientos <- (paste(matriz[initialrow:nrow(matriz),treatment1col],matriz[initialrow:nrow(matriz),treatment2col], sep="&"))
  vectortratamientos<-factor(vectortratamientos, levels = unique(vectortratamientos))
  #### Consistentes ####
  datos <- split(as.data.frame(matriz[initialrow:nrow(matriz),initialcolumn:ncol(matriz)]),vectortratamientos,drop=T) 
  maxsignificantvalues <-lapply(datos,function(x) colSums(x !=0)/nrow(x)) 
  consist <- do.call(pmax,maxsignificantvalues)
  consist <- consist ==1 
  #### Cuantitativo ####
  minreplicates <- round(nrow(matriz[initialrow:nrow(matriz),])*threshold,0) 
  cuantit <- colSums(matriz[initialrow:nrow(matriz),initialcolumn:ncol(matriz)] !=0) 
  cuantit <- cuantit >=minreplicates 
  #### Seleccion ####
  seleccion <- consist|cuantit 
  cat("\nVariable Selection based on consistency. Variable was present at least on",minreplicates, " cases, or all replicates of a treatment" )
  cat("\nInitial Variables: ", ncol(matriz[,initialcolumn:ncol(matriz)]), " Selected Variables: ", length(seleccion[seleccion==T]), " Removed Variables: ", length(seleccion[seleccion==F]))
  matrizfiltrada <- cbind(as.data.frame(matriz[,1:initialcolumn-1]),as.data.frame(matriz[,initialcolumn:ncol(matriz)][,seleccion==T])) 
  colnames(matrizfiltrada)[1:initialcolumn-1]<-colnames(matriz)[1:initialcolumn-1]
  return(matrizfiltrada)
  
}
