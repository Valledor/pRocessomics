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


RemoveEmptyColumns_import <-function(matriz, initialrow, initialcolumn, datasetname){
  matrizfiltrada <-matriz[initialrow:nrow(matriz),initialcolumn:ncol(matriz)]
  colnamesorig <-colnames(matrizfiltrada)
  vectorseleccion <- colSums(matrizfiltrada,na.rm=TRUE)==0
  if(length(vectorseleccion)==0){matrizfiltrada<-matrizfiltrada}
  else {matrizfiltrada <- matrizfiltrada[,vectorseleccion==F]
  if(length(vectorseleccion[which(vectorseleccion)>0]))
  {
    cat(paste("\n\n WARNING: In", datasetname, "dataset, some columns were empty and removed for further analyses:\n\n", paste(colnamesorig[vectorseleccion==T],collapse=", ")))    
    cat("\n\nHint: Sometimes you can get datasets with empty columns, but this is not usual. If you don't feel that this is normal, please, re-check your data.\n\n")
  }
  }
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
