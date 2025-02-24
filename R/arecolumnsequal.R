#
# Copyright 2019 Luis Valledor / Laura Lamelas
# Last revision 31.05.2024
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

arecolumnsequal <-function(datasetlist,column){
  samplenames <- lapply(datasetlist, function(x) x[,column])
  samplenames <- do.call(cbind,samplenames)
  samplenamesequal <- c()
  for (i in 1:length(datasetlist)){
    aux <- identical(samplenames[,1],samplenames[,i])
    samplenamesequal <- c(samplenamesequal, aux)
  }
  return(samplenamesequal)
}


