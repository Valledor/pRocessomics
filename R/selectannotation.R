#
# Copyright 2019 Luis Valledor / Laura Lamelas
# Last revision 30.05.2019
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

#### Select annotation file
selectannotation <- function(annotationfilename=NULL){
  anotlist <-c("foo")
  if(is.null(annotationfilename)==FALSE){
    datalist <-get(annotationfilename,envir=.GlobalEnv)  
    cat("\nYou have selected",annotationfilename, "annotation file.")
    writeLines(annotationfilename, con = stdout(), sep = "\n", useBytes = FALSE) 
    
  }
  if(is.null(annotationfilename)==TRUE){
    anotlists <- Filter(function(x) "pRoAnnot" %in% class(get(x)), ls(envir = .GlobalEnv))
    if(length(anotlists)==0) stop("No valid pRoAnnot objects were found. Please re-import your annotations using importannotation() function.")
    texts_wizard("We have detected the following annotation files available. Please enter the number of the annotation file do you want to use and press [enter].")
    annotationfilename <- utils::select.list(anotlists, preselect = NULL, multiple = F, title = NULL, graphics = FALSE)
    writeLines(annotationfilename, con = stdout(), sep = "\n", useBytes = FALSE) 
    anotlist <-get(annotationfilename,envir=.GlobalEnv)
  }
  #Now we check if we have a list (pRoAnnot object) or a list of lists, in this case we ask which list 
  if(class(anotlist[[1]]) %in% c("pRoAnnot","list")==TRUE){
    cat("\nThis annotation_file contains multiple choices. Please select which transformation do you want to employ")    
    trflist<-names(anotlist)
    mylist<- utils::select.list(trflist, preselect = NULL, multiple = F, title = NULL, graphics = FALSE)
    writeLines(mylist, con = stdout(), sep = "\n", useBytes = FALSE) 
    anotlist <-anotlist[[mylist]]
  }
  tempannotationfilename <- list(tempannotationfilename=annotationfilename)
  list2env(tempannotationfilename,envir = .GlobalEnv)
  return(anotlist)
}
