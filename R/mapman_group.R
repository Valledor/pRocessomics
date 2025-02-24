#
# Copyright 2019 Luis Valledor / Laura Lamelas
# Last revision 01.2021
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

#' @name mapman_group
#' @title Mapman group
#' @description A function to sort and group the variables according to Mapman categories or custom annotations
#' @usage mapman_group(datalist, annotation, initialrow = 1, initialcolumn = 3, 
#' treatment1col = 1, treatment2col = 2, treatment = 1,
#' group = TRUE, omiclevel = NULL, threshold = 1, stats = FALSE, scalevalues = FALSE)
#' @param datalist List with different preprocessed omic levels. POL class object.
#' @param annotation "annot" class object containing variables descriptions according to mercator mapman, generated from pRocessomics::importannotation() function
#' @param initialrow First row of numeric data within datasets.
#' @param initialcolumn First column of numeric data within datasets.
#' @param treatment1col Column number in which treatment 1 is indicated.
#' @param treatment2col Column number in which treatment 2 is indicated (if applicable).
#' If there is only one treatment indicate the same column as treatment1col.
#' @param treatment numeric. Set 1 to split the dataset according to treatment1, 2 to split the
#' dataset according to treatment 2 and 3 to split the dataset according to the combination of both treatments.
#' @param group Boolean to indicating whether replicates from same treatment should be grouped. If set to positive average values will be shown
#' @param omiclevel Character vector indicating the quoted name or names of the omiclevel(s) to be analyzed
#' @param threshold Number indicating number of minimum variables in a Mapman category to consider that category. Default 1.
#' @param stats Boolean indicating whether stats calculations as standard deviation and pvalues should be returned with the mapman table
#' @param scalevalues Boolean indicating whether the values of the variables should be scaled to carry out the pvalue calculations. It will be only considered if stats = TRUE
#' @return Mapman grouped table with sums of abundances according to Mapman categories
#' @author Luis Valledor and Laura Lamelas

#' @importFrom stats sd aov anova
#' @export

mapman_group <- function(datalist, annotation, initialrow = 1, initialcolumn = 3, treatment1col = 1, treatment2col = 2, treatment = 1,
                             group = TRUE, omiclevel = NULL, threshold = 1, stats = FALSE, scalevalues = FALSE) {
  options(warn=-1)
  #### Initial checks ####
  if(!inherits(datalist, "POL")) stop("A POL object is required. Please run pre-processing steps")
  if(is.null(annotation)) stop("\nAnnotation matrix is required. Please run MapMan to annotate your variables\nVisit https://mapman.gabipd.org for more information.")
  datasetnames <- names(datalist)
  if(is.null(omiclevel)) omiclevel <- datasetnames
  if(any(is.na(omiclevel))) omiclevel <- datasetnames
  sapply(omiclevel, function(x) if(x %in% datasetnames == FALSE) stop(paste("Error: Please select a valid datasetname (",paste(datasetnames,collapse = " ")),")",sep=""))
  texts_wizard("\n\nGROUPING VARIABLES ACCORDING TO MAPMAN\n\nSingle thread computing. It may take a while...\n")
  
  #### Selection of omic levels and extraction of meaningful rows and cols ####
  #Remove non significant rows and cols    
  filtdatalist <- lapply(datalist, function(x) x[initialrow:nrow(x),initialcolumn:ncol(x)])    
  #we generate a data matrix with selected datasets
  if(length(omiclevel)>1){  
    variablesformapman <-do.call(cbind,filtdatalist[omiclevel])
    variablenames <- lapply(filtdatalist[omiclevel], function(x) colnames(x))
    variablenames <-unlist(variablenames)
    colnames(variablesformapman) <- variablenames
    variablesformapman<-RemoveEmptyColumns(variablesformapman,1,1) #This avoids strange behaviors
  }
  if(length(omiclevel)==1){
    variablesformapman <-filtdatalist[[omiclevel]]
    variablenames <-colnames(variablesformapman)
    variablesformapman<-RemoveEmptyColumns(variablesformapman,1,1) #This avoids strange behaviors
  }
  
  #### Preparation of annotations and groups ####
  IDENTIFIER <- annotation$IDENTIFIER
  annotation<-subset(annotation, IDENTIFIER %in% variablenames, select=colnames(annotation))
  bincodeaux <- seq(from = 3, to = ncol(annotation), by = 2)
  nameaux <- seq(from = 4, to = ncol(annotation), by =2)
  
  # Check annotations.
  if(grepl(".", annotation[1,3], fixed = TRUE)==T){ #that should be impossible but I feel like keeping it anyway
   for(i in bincodeaux){
     annotation[, i] <- as.numeric(do.call(rbind, lapply(strsplit(as.vector(annotation[, i]), "\\."), function(x) x[[1]])))
   }
   
   for(i in nameaux){
     annotation[, i] <- as.character(do.call(rbind, lapply(strsplit(as.vector(annotation[, i]), "\\."), function(x) x[[1]])))
   }
  }
  
  #### reshaping of annotation matrix ####
  anf <-c()
  for(i in nameaux){ #unfolding the annot.matrix
    anf2<-annotation[,c(1,i)][!is.na(annotation[,i]),]
    colnames(anf2) <- c("IDENTIFIER","NAME")
    anf<- rbind(anf,anf2)
    }
  singles <- which(anf[,2]=="")
  if(length(singles)!=0) anf <- anf[-singles,]
  
  #### mapman groups present in annotations ####
  groupsmm <-unique(anf[,2])
  
  #### sums per groupsmm ####
  groupedmatrix <- c()
  elemxbin <- c()
  for(i in groupsmm){
    misvar <-as.vector(anf[anf[,2]==i,1])
    elemxbin<-c(elemxbin,length(misvar))
    if(length(misvar)<threshold) next
    if(length(misvar)==1)groupedmatrix <- cbind(groupedmatrix,variablesformapman[,colnames(variablesformapman)==misvar])
    if(length(misvar)>1)groupedmatrix <- cbind(groupedmatrix,rowSums(variablesformapman[,colnames(variablesformapman)%in%misvar]))
  }
  colnames(groupedmatrix)<-groupsmm
  names(elemxbin) <- groupsmm
  
  #### define treatments ####
  splitvector <- as.vector(buildtreatments(datalist,c(initialrow, initialcolumn, treatment1col, treatment2col))[,treatment])
  
  #### Export results when not grouping #### 
  if (group == F) {
    result <- list(t(groupedmatrix), elemxbin, splitvector)
    names(result) <- c("mean","number_bin_elements", "treatments")
    class(result) <- "MMO"
    texts_wizard("\nmapman_group: done!")
    return(result)
  }
  
  #### Grouping by treatment ####
  if(group == T) {
    tsv <- unique(splitvector)
    ttgm <- c()
    ttgmsd <- c()
    for(i in tsv){
      ttgm <- cbind(ttgm,colMeans(groupedmatrix[which(splitvector==i),]))
      ttgmsd <- cbind(ttgmsd,sqrt(diag(var(groupedmatrix[which(splitvector==i),]))))
    }
    colnames(ttgm)<-tsv
    colnames(ttgmsd)<-tsv
    pvalues<-NA
    
    if(stats==T){
      
      if(scalevalues == T) groupedmatrix <- scale(groupedmatrix)
      groupedmatrix<-as.data.frame(groupedmatrix)
      for (i in 1:ncol(groupedmatrix)) {
        res.anova = stats::aov(groupedmatrix[, i] ~ (splitvector), data = groupedmatrix,na.rm = TRUE)
        pvalues <- c(pvalues, stats::anova(res.anova)$'Pr(>F)'[1:1])
        names(pvalues)[1:length(pvalues)] <- colnames(groupedmatrix)
        
      }
    }
    
    result <- list(ttgm, ttgmsd, elemxbin, pvalues, splitvector)
    names(result) <- c("mean", "sd", "number_bin_elements", "pvalues", "treatments")
    class(result) <- "MMOG"
    texts_wizard("\nmapman_group: done!")
    return(result)
  }
}


