#
# Copyright 2019 Luis Valledor / Laura Lamelas
# Last revision 04.11.2019
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
RemoveEmptyColumns <-function(matriz, initialrow, initialcolumn){
  matrizfiltrada <-matriz[initialrow:nrow(matriz),initialcolumn:ncol(matriz)]
  initial<-ncol(matrizfiltrada)
  vectorseleccion <-which(colSums(matrizfiltrada,na.rm=TRUE)==0)
  end<-length(vectorseleccion)
  cat(paste("\nOut of ",initial,"initial variables, ",end, "have been found to be empty and dropped from the analysis"))
  if(length(vectorseleccion)==0){
    return(matriz)}
  else {matrizfiltrada <- matrizfiltrada[,-vectorseleccion]}
  if(initialrow==1){
    matrizfiltrada <- cbind(matriz[,1:initialcolumn-1],matrizfiltrada)
    return(matrizfiltrada)
  } else if(initialrow>=2) {
    matrizfiltrada <- cbind(matriz[,1:initialcolumn-1],matrizfiltrada)
    matrizfiltrada <- rbind(matriz[1:initialrow-1,],matrizfiltrada)
    
    return(matrizfiltrada)
  } else {
    return(cat("\nError. Please revise your data input."))
  }
}
