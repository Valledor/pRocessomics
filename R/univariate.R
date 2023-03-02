#
# Copyright 2019 Luis Valledor / Laura Lamelas
# Last revision 01.2023
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

#' @name univariate
#' @title Univariate statistical analysis.
#' @description This function performs univariate statistical analyses such as 
#' parametrics as t test or ANOVA and non-parametric tests as Wilcox test and 
#' Kruskal Wallis, there is also the posibility to adjust obtained pvalues, 
#' FDR and post hoc analyses.
#' @usage univariate(datalist,initialrow=1,initialcolumn=3, treatment1col=1,
#' treatment2col=2, treatment=1, parametric=TRUE, posthoc=TRUE, FDR=TRUE, 
#' round=6, annotatefile=NULL)
#' @param datalist List with different preprocessed omic levels. POL class object.
#' @param initialrow First row of numeric data within datasets.
#' @param initialcolumn First column of numeric data within datasets.
#' @param treatment1col Column number in which treatment 1 is indicated.
#' @param treatment2col Column number in which treatment 2 is indicated (if applicable).
#' If there is only one treatment indicate the same column as treatment1col.
#' @param treatment numeric. Set 1 to split the dataset according to treatment1, 2 to split the
#' dataset according to treatment 2 and 3 to split the dataset according to the combination of both treatments.
#' @param parametric logical, this argument indicates whether the employed test 
#' should be or not parametric, depending on the number of available treatments if parametric a t test or ANOVA 
#' test will be employed, for non parametric data, Wilcox test or kruskal Wallis test will be employed.
#' @param posthoc logical, this argument indicates whether to do a posthoc test, if parametric=TRUE a 
#' TukeyHSD test will be done, for non-parametric analysis Dunn test.
#' @param FDR logical, if TRUE adjusted pvalues will be provided according to Bonferroni method
#' @param round numeric, number of desired decimal digits in the output
#' @param annotatefile optional, pRoAnnot class object created with pRocessomics::importannotation function
#' @details The objective of univariate function is to perform statistical analyses
#' in order to check variability across treatments.
#' This function requires a POL object, previously defined by using
#' preprocess_omic_list or preprocess_wizard and, if desired, transformdata/transformation_wizard or featureselection/featureselection_wizard. 
#' Two different approaches are implemented within this function, the first related 
#' to parametriciy, and the second taking into account the number of treatments to 
#' compare. 
#' @return  List containing mean, standard deviation and selected stats from each omic layer, 
#' UNIVAR class object.
#' @author Luis Valledor and Laura Lamelas
#' @seealso transformdata and featureselection
#' 
#' @export



univariate <-function(datalist, initialrow=1, initialcolumn=3, treatment1col=1, 
                      treatment2col=2, treatment=1,parametric=TRUE, 
                      posthoc=TRUE, FDR=TRUE,round=6,annotatefile=NULL){
  options(warn=-1)
  
  #suppressPackageStartupMessages(require(dplyr))
  #suppressPackageStartupMessages(require(plyr))
  #suppressPackageStartupMessages(require(FSA))
  
  # Initial checks
  if(!inherits(datalist, "POL")) stop("A POL object is required. Please run pre-processing steps")
  datasetnames<-names(datalist)
  
  #Empty column removal
  datalist <- lapply(datalist, function(x) RemoveEmptyColumns(x, initialrow=initialrow,initialcolumn=initialcolumn))

  #numeric values
  listadatos2 <- lapply(datalist, function(x) x[initialrow:nrow(x),initialcolumn:ncol(x)])

  #grouping by treatment
  if(treatment==1) splitvector <- as.vector(datalist[[1]][,treatment1col])
  if(treatment==2) splitvector <- as.vector(datalist[[1]][,treatment2col])
  if(treatment==3) splitvector <- as.vector(paste(datalist[[1]][,treatment1col], datalist[[1]][,treatment2col],sep="&"))
  
  stattest <- c()
  if(length(unique(splitvector))==1) stop("Two or more treatments are required")
  if(length(unique(splitvector))==2) stattest<-"t_or_u"
  if(length(unique(splitvector))>=3) stattest<-"anova_or_kruskal"
  
  splitvector<-factor(splitvector,levels=unique(splitvector))
  
  if(stattest=="t_or_u"){
    if(parametric==TRUE){
      pvaluelist<-lapply(datalist, function(x) testT(x,splitvector,initialcolumn,initialrow))  
      #pvaluelist<-lapply(pvaluelist, function(x) apply(x,2,function(y) format(y,scientific = FALSE)))
      if(FDR==TRUE) pvaluelist <- lapply(pvaluelist, function(x) qvalues(x, round))
    }
    if(parametric==FALSE){
      pvaluelist<-lapply(datalist, function(x) testU(x,splitvector,initialcolumn,initialrow))  
      #pvaluelist<-lapply(pvaluelist, function(x) apply(x,2,function(y) format(y,scientific = FALSE,digits = round)))
      if(FDR==TRUE) pvaluelist <- lapply(pvaluelist, function(x) qvalues(x, round))
    }
  
  }
  
  if(stattest=="anova_or_kruskal"){
    if(parametric==TRUE){
      pvaluelist<- lapply(listadatos2,function(y) apply(y,2,function(x) funcionanova(x,splitvector,posthoc,round)))
      #pvaluelist<-lapply(pvaluelist, function(x) apply(x,2,function(y) format(y,scientific = FALSE)))
      if(FDR==TRUE) {pvaluelist <- lapply(pvaluelist, function(x) qvalues(x, round))
      }
    }
    if(parametric==FALSE) {
        pvaluelist<- lapply(listadatos2,function(y) apply(y,2,function(x) kruskalwallis(x,splitvector,posthoc,round)))
        #pvaluelist<-lapply(pvaluelist, function(x) apply(x,2,function(y) format(y,scientific = FALSE)))
        if(FDR==TRUE) {pvaluelist <- lapply(pvaluelist, function(x) qvalues(x, round))
        }
    }
  }
  
  #tablas de medias y SDs
  meansdlist <- lapply(listadatos2,function(x) meansd(x,splitvector,round))

  #componemos tabla final
  compositelist <- meansdlist
 
  for(i in 1:length(compositelist)){
    if (FDR==F&posthoc==F) {
      #pvaluelist[[i]]<-apply(pvaluelist[[i]],c(1,2), function(x) format(x,scientific = FALSE))
      pvaluelist[[i]]<-as.data.frame(pvaluelist[[i]])
      colnames(pvaluelist[[i]])<-"p-value"
      #pvaluelist[[i]]<-pvaluelist[[i]][order(pvaluelist[[i]]$`p.value`),]#comentado porque me quita los rownames
      compositelist[[i]] <- cbind(compositelist[[i]], pvaluelist[[i]])}
    else {compositelist[[i]] <- cbind(compositelist[[i]], t(pvaluelist[[i]]))}
  }
  #aqui esta composite list terminada
  if(!is.null(annotatefile)){
    compositelist<-lapply(compositelist,function (x) do.annotate(annotatefile,x))
  }
  results <- list(meansd=meansdlist, pvalue=pvaluelist, composite=compositelist)
  class(results) <- "UNIVAR"
  #univar_plot(results,annotatefile)
 return(results)
}

