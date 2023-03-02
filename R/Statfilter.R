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

#' @importFrom stats anova aov kruskal.test p.adjust

Statfilter <- function(matriz, initialrow, initialcolumn, treatment1col, treatment2col, treatment, threshold=0.25, parametric=TRUE, stat="p",abovebelow="Below"){
  datos<-matriz[initialrow:nrow(matriz),initialcolumn:ncol(matriz)] 
  if(treatment==1) tratamientos <- (matriz[initialrow:nrow(matriz),treatment1col])
  if(treatment==2) tratamientos <- (matriz[initialrow:nrow(matriz),treatment2col])
  if(treatment==3) tratamientos <- (paste(matriz[initialrow:nrow(matriz),treatment1col],matriz[initialrow:nrow(matriz),treatment2col], sep="&"))
  tratamientos<-factor(tratamientos, levels = unique(tratamientos))
  if(parametric==TRUE){
    pvalores<-sapply(datos, function(x) stats::anova(stats::aov(x~tratamientos))$'Pr(>F)'[1])
    message<-"ANOVA"
  } else {
    pvalores<-sapply(datos, function(x) stats::kruskal.test(x~tratamientos)$p.value)
    message<-"Kruskal"
  }
  if(stat=="p"){
    if(abovebelow=="Above") {
      seleccion <- pvalores>=threshold
    }
    if(abovebelow=="Below"){
      seleccion <- pvalores<=threshold
    }

  } else if(stat=="q") {
    qvalores <-stats::p.adjust(pvalores, method = "BH")
    
    if(abovebelow=="Above") {
      seleccion <- qvalores>=threshold
    }
    if(abovebelow=="Below"){
      seleccion <- qvalores<=threshold
    }
  } else {
    return(cat("\nError. stat = c('p','q'), Please select mode p or q value"))
  }
  matrizfiltrada<-as.data.frame(datos[,seleccion])
  colnames(matrizfiltrada)<-colnames(datos[seleccion])

    dog <-paste("\nStat filtering: Out of",length(seleccion),"variables,",length(seleccion[seleccion==T]),"showed a", stat, "value", abovebelow, threshold,"threshold.",sep=" ")
    texts_wizard(dog)
    
    
  if(initialrow==1){
    matrizfiltrada <- cbind(matriz[,1:initialcolumn-1],matrizfiltrada)
    colnames(matrizfiltrada)[1:initialcolumn-1]<-colnames(matriz)[1:initialcolumn-1]
    return(matrizfiltrada)
  } else if(initialrow>=2) {
    matrizfiltrada <- cbind(matriz[,1:initialcolumn-1],matrizfiltrada)
    colnames(matrizfiltrada)[1:initialcolumn-1]<-colnames(matriz)[1:initialcolumn-1]
    matrizfiltrada <- rbind(matriz[1:initialrow-1,],matrizfiltrada)
    return(matrizfiltrada)
  } else {
    return(cat("\nError. Please revise your data input."))
  }
}

