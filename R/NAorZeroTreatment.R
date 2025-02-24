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

NAorZeroTreatment <- function(submatriz,threshold,initialcolumn){
  submatrix <- submatriz[,initialcolumn:ncol(submatriz)]
  replicates <- nrow(submatrix)
  NACount <- apply(submatrix,2, function(x) sum(is.na(x)))
  Isa0 <- NACount <= (replicates*threshold)
  for(i in 1:length(Isa0)){
    if(Isa0[i]==FALSE){submatrix[,i][is.na(submatrix[,i])]<-0}
  }
  submatriz[,initialcolumn:ncol(submatriz)] <- submatrix
  return(submatriz)
}
