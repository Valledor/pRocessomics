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

#' @importFrom stats aov anova TukeyHSD


funcionanova <- function(vector,tratamiento,posthoc=TRUE,round=5){
  res.anova=stats::aov(vector~as.factor(tratamiento), na.rm=TRUE)
  pvalor <- stats::anova(res.anova)$'Pr(>F)'[1:1]
  names(pvalor)[1] <- "p-value"
  valordevuelta <-(pvalor)
  #valordevuelta<-format(valordevuelta,scientific = FALSE)
  
  if(posthoc==TRUE){
    res.tukey <- stats::TukeyHSD(res.anova)
    pvalores.tukey <-res.tukey[[1]][,4]
    valordevuelta<-c(pvalor,pvalores.tukey)
    valordevuelta <- (round(valordevuelta,round))
    return(valordevuelta)
  }
  
  valordevuelta<-(as.numeric(as.character(valordevuelta)))
  valordevuelta <- t(round(valordevuelta,round))
  
  return(valordevuelta)
}
