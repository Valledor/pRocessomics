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

#' @importFrom stats p.adjust

qvalues <- function(x,round=5){
  # This transform vector in data frame, when only p values and not Kruskal are requested
  if(is.vector(x)){
    x <-t(as.data.frame(x))
    row.names(x) <- "p-value"
  }
  if(ncol(x)==1){
    x <- t(x)
  }
  qval<-stats::p.adjust(x[1,],method = "BH")
  qval<-round(qval,round)
  #print(nrow(x))
  if(nrow(x)>=2) {y <- rbind(x[1,],qval,x[2:nrow(x),])}
  if(nrow(x)==1) {y <- rbind(x[1,],qval)}
  
  row.names(y)[1:2] <- c(rownames(x)[1],"q-value")
  return(y)
}
