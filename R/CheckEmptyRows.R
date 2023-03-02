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


CheckEmptyRows <-function(matriz, initialrow, initialcolumn){
  matrizfiltrada <-matriz[initialrow:nrow(matriz),initialcolumn:ncol(matriz)]
  colnamesorig <-colnames(matrizfiltrada)
  vectorseleccion <- rowSums(matrizfiltrada,na.rm=TRUE)==0
  return(vectorseleccion)
}
