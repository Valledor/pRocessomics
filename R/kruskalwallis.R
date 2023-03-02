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

#' @importFrom stats kruskal.test
#' @importFrom FSA dunnTest

kruskalwallis <- function(vector,tratamiento,posthoc=TRUE,round=5){
  reskruskal<-stats::kruskal.test(vector~as.factor(tratamiento))
  pvaloresKruskal<-reskruskal$p.value
  names(pvaloresKruskal)[1] <- "p-value"
  valordevuelta <-as.matrix(pvaloresKruskal)
  if(posthoc==TRUE){
    #require("FSA")
    dt<-FSA::dunnTest(vector~as.factor(tratamiento),method = "bonferroni")
    pvaldunn<-dt[[2]][,3]
    names(pvaldunn)<-dt[[2]][,1]
    valordevuelta<-c(pvaloresKruskal,pvaldunn)
    valordevuelta <- round(valordevuelta,round)
    return(valordevuelta)
  }
  
  valordevuelta <- t(round(valordevuelta,round))
  return(valordevuelta)
}
