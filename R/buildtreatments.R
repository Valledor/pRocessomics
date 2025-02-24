#
# Copyright 2019 Luis Valledor / Laura Lamelas
# Last revision 26.03.2019
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
buildtreatments <- function(datalist,settings){
  if(is.na(settings[4])==TRUE){
    uniquetreatment <- as.vector(datalist[[1]][settings[1]:nrow(datalist[[1]]),settings[3]])
    treatments <- cbind(uniquetreatment, uniquetreatment, uniquetreatment)
  } else {
    treatments <- cbind(as.vector(datalist[[1]][settings[1]:nrow(datalist[[1]]),settings[3]]), as.vector(datalist[[1]][settings[1]:nrow(datalist[[1]]),settings[4]]), paste(as.vector(datalist[[1]][settings[1]:nrow(datalist[[1]]),settings[3]]), as.vector(datalist[[1]][settings[1]:nrow(datalist[[1]]),settings[4]]),sep="&"))
  }
  treatments <- apply(treatments,2,function(x) { factor(x, levels=unique(x))
    return(x)})
  treatments<-as.data.frame(treatments)
return(treatments)
  }
