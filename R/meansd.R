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

meansd <- function(table,splitvector,decimals=3){
  listpertreatments <- split(table,splitvector,drop=T)
  means<-lapply(listpertreatments, colMeans)
  means<-as.data.frame(do.call(cbind, means))
  sds <- lapply(listpertreatments, function(x) apply(x,2,stats::sd))
  sds<-as.data.frame(do.call(cbind,sds))
  means<-round(means,decimals)
  sds<-round(sds,decimals)
  plusminussign<- " \u00B1 "
  Encoding(plusminussign)<-"UTF-8"
  plusminus <- rep(c(plusminussign),nrow(means))
  outputtable<-c()
  for(i in 1:ncol(means)){
    temptable<-cbind.data.frame(means[,i],plusminus,sds[,i])
    colnames(temptable)<- c(paste("Mean-",(colnames(means)[i]),sep="")," ",paste("SD-",(colnames(means)[i]),sep=""))
    outputtable <-c(outputtable,temptable)
    outputtable <- as.data.frame(outputtable)
  }
  row.names(outputtable)<-colnames(table)
  

  return(outputtable)
}
