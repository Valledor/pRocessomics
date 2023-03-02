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

#' @importFrom plyr count
#' @importFrom stats wilcox.test

testU<-function(matriz,splitvector, initialcolumn, initialrow){
  #replicas por tratamiento?
  nrep<-plyr::count(splitvector)[1,2]
  #haria falta un check de que todos los tratamientos de todas las variables tienen el mismo numero de replicas
  #Loop
  p.fin.t<-c()
  for(j in initialcolumn:ncol(matriz)){
    a<-stats::wilcox.test(matriz[initialrow:(initialrow+nrep-1),j],matriz[initialrow+nrep:nrow(matriz),j])$p.value
    p.fin.t<-c(p.fin.t,a)
    j=j+1
  }
  p.fin.t<-as.data.frame(t(p.fin.t))
  colnames(p.fin.t)<-colnames(matriz[initialcolumn:ncol(matriz)])
  rownames(p.fin.t)<-"p-value"
  return(t(p.fin.t))
}
