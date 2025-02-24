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

#' @name featureselection
#' @title Feature Selection in Omic Datasets
#' @description This function processes a POL object to be analyzed in order to select featured variables
#' @usage featureselection(datalist, initialrow, initialcolumn, 
#' treatment1col, treatment2col, treatment, method = NULL, 
#' threshold = NULL, parametric = TRUE, stat = "p", abovebelow = "Below")
#' @param datalist List with different preprocessed omic levels. pRoDS class object.
#' @param initialrow First row of numeric data within datasets.
#' @param initialcolumn First column of numeric data within datasets.
#' @param treatment1col Column number in which treatment 1 is indicated.
#' @param treatment2col Column number in which treatment 2 is indicated (if applicable).
#' If there is only one treatment indicate the same column as treatment1col.
#' @param treatment numeric. Set 1 to split the dataset according to treatment1, 2 to split the
#' dataset according to treatment 2 and 3 to split the dataset according to the combination of both treatments.
#' @param method Method for selecting variables: "IQR"- interquartile range; "Stat" - ANOVA/Kruskal Wallis; "VarCoef" - Coeficient of variation.
#' @param threshold Desired threshold to apply to the data
#' For "IQR" : % over the average IQR; For "Stat": p or q value; For "VarCoef" : Coeficient of variation
#' @param parametric boleean (TRUE, FALSE). If TRUE ANOVA test will be used, if FALSE Kruskal Wallis.
#' @param stat Indicates which statistic will be used p-values,"p", or q-values, "q". Only required when method = "Stat".
#' @param abovebelow character, set to "Above" for keep those variable- values above the threshold or "Below" to keep those variables below the threshold
#' @details The objective of \code{featureselection} is removing those variables who do not greatly change between treatments, known to later introduce noise in multivariate analyses reducing the power of these approaches. This function requires a \code{POL} object, previously defined by using \code{preprocess_omic_list} and, if desired, \code{transformdata}. Two different approaches are implemented within this function: first one is based on selecting those variables that have intercuartilic ranges, IQR, greater than average within all studied variables. In this case the \code{threshold} indicates the \% above average IQR and only variables with greater or equal IQR will be selected. The second strategy is based on applying ANOVA or Kruskal Wallis and selected those variables exhibiting a p or q value smaller or equal than \code{threshold}.
#' @return returns a "POL" object with the filtered data
#' @author Luis Valledor and Laura Lamelas
#' @export
#' 

featureselection<-function(datalist, initialrow, initialcolumn, treatment1col, treatment2col, treatment, method=NULL, threshold=NULL, parametric=TRUE,stat="p",abovebelow="Below"){
  if(!inherits(datalist, "POL")) stop("A POL object is required. Please run pre-processing steps")
  
  datasetnames <- names(datalist)
  if(!is.null(method)){
    texts_wizard("\nFEATURE SELECTION")
    texts_wizard(paste("\n",method, " filtering: Single processor core will be used. It may take a while...",sep=""))
    datalist<-mapply(function(x,y,z,a,b,c){
      feature_selection_list_helper(x,initialrow=initialrow,initialcolumn=initialcolumn,treatment1col=treatment1col,treatment2col=treatment2col,
                                    treatment=treatment,method=y,threshold=z,parametric=a,stat=b,datasetname=c,abovebelow=abovebelow)
    },datalist,method,threshold,parametric,stat,datasetnames,SIMPLIFY = FALSE)
  }
  class(datalist)<-"POL"
  return(datalist)
}
