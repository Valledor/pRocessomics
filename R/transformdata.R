#
# Copyright 2019 Luis Valledor / Laura Lamelas
# Last revision 01.2023
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
#' @name transformdata
#' @title Data Transformation
#' @description A function to transform omic data
#' @usage transformdata(datalist, initialrow=1, initialcolumn=3, treatment1col = 1, 
#' treatment2col = 2, treatment = 3, transf=NULL)
#' @param datalist List with different preprocessed omic levels. pRoDS class object.
#' @param initialrow First row of numeric data within datasets.
#' @param initialcolumn First column of numeric data within datasets.
#' @param treatment1col Column number in which treatment 1 is indicated.
#' @param treatment2col Column number in which treatment 2 is indicated (if applicable).
#' If there is only one treatment indicate the same column as treatment1col.
#' @param treatment numeric. Set 1 to split the dataset according to treatment1, 2 to split the
#' dataset according to treatment 2 and 3 to split the dataset according to the combination of both treatments.
#' @param transf Character vector, where transf and datalist should have the same length determining the method used to transform each omiclevel. If only one value is provided as transf="Log" all the omiclevel will be log transformed, available options are:
#' \itemize{
#' \item z <- data scaling
#' \item z-centered <- data centering and scaling
#' \item Log10 <- logarithmic transformation *please note negative values will yield a NaN
#' \item Log10plus1 <- log(data+1) *please note negative values will yield a NaN
#' \item box-cox-lambda <- box cox lambda transformation
#' \item sqrt <- square root *please note negative values are not allowed
#' \item cubic <- cubic root
#' \item arcsine <- arcsine *please note maximum absolute value is 100
#' }
#' @return A list containing the transformed dataset
#' @author Luis Valledor and Laura Lamelas
#' @seealso transformation_wizard

#' 
#' @export

transformdata <- function(datalist, initialrow=1, initialcolumn=3, treatment1col = 1, treatment2col = 2, treatment = 3, transf=NULL){
  #### Initial checks ####
  if(!inherits(datalist, "POL")) stop("A POL object is required. Please run pre-processing steps")
  
  # Set and check transformation 
  if(!is.null(transf)) if(length(transf)!=1&length(transf)!=length(datalist)) stop("A vector containing one-common- or n elements (where n is the number of datasets) indicating transformation method(s) is expected")
  if(!is.null(transf)) if(FALSE %in% (transf %in% c('z','z-centered','Log10', 'Log10plus1', 'sqrt', 'cubic', 'arcsine', 'box-cox-lambda','none'))==TRUE) stop("Please select a valid transformation method")
  if(!is.null(transf)) if(length(transf)==1) transf <- rep(transf, length(datalist))
  
  
  datasetnames <-names(datalist)
  
  #### Transformation ####
  if(!is.null(transf)){
    texts_wizard("\nTRANSFORMATION OF DATASETS\n")
    texts_wizard("Single processor core will be used. It may take a while...\n")
    cat(paste(transf, " transformation method will be used for ", names(datalist)," dataset\n",sep=""))
    datalist<-mapply(function(x,y,z){
      Trnsfrm(x,initialrow=initialrow,initialcolumn=initialcolumn,treatment1col = treatment1col, treatment2col = treatment2col, treatment = treatment,datasetname = z,transformation = y)
    },datalist,transf,datasetnames,SIMPLIFY = FALSE)
  }
  cat("\ndone!")
  
  class(datalist)<-"POL"
  
  return(datalist)
}
