#
# Copyright 2019 Luis Valledor / Laura Lamelas
# Last revision 25.03.2019
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
#' @importFrom utils select.list

#### Select datalist
selectdataset <- function(datasetname=NULL){
  datalist <-c("foo")
  if(is.null(datasetname)==FALSE){
    datalist <-get(datasetname,envir=.GlobalEnv)  
  dog<-as.character("\nYou have selected",datasetname, "dataset.")
  texts_wizard(dog)
    }
  if(is.null(datasetname)==TRUE){
  datalists <- Filter(function(x) "POL" %in% class(get(x)), ls(envir = .GlobalEnv))
  if(length(datalists)==0) stop("No valid POL objects were found. Please pre-process your data first.")
  texts_wizard("We have detected the following datasets available. Please enter the number of the dataset do you want to process and press [enter].")
  datasetname <- utils::select.list(datalists, preselect = NULL, multiple = F, title = NULL, graphics = FALSE)
  datalist <-get(datasetname,envir=.GlobalEnv)
  }
  #Now we check if we have a list (POL object) or a list of lists, in this case we ask which list 
  if(class(datalist[[1]]) %in% c("POL","list")==TRUE){
    texts_wizard("\nThis dataset contains multiple choices. Please select which transformation do you want to employ")    
    trflist<-names(datalist)
    mylist<- utils::select.list(trflist, preselect = NULL, multiple = F, title = NULL, graphics = FALSE)
    datalist <-datalist[[mylist]]
  }
  tempdatasetname <- list(tempdatasetname=datasetname)
  list2env(tempdatasetname,envir = .GlobalEnv)
  return(datalist)
  }
