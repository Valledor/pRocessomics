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

#' @importFrom  impute impute.knn
#' @importFrom doParallel registerDoParallel

NAorZero_Imput <- function(matriz, initialrow, initialcolumn, treatment1col, treatment2col, treatment=1, threshold=0.25, imputation="RF",k=3, cores=2,datasetname=NULL){
  cat(paste("\nProcessing ",datasetname," dataset",sep=""))
  #suppressMessages(suppressWarnings(require("plyr")))
  datos<-matriz[initialrow:nrow(matriz),initialcolumn:ncol(matriz)]
  initial0s <- sum(datos==0,na.rm=TRUE)
  initialNAs <-sum(is.na(datos))
  variablesinitial0s <- length(apply(apply(datos,2,function(x) x==0),2,any,na.rm=T)[apply(apply(datos,2,function(x) x==0),2,any,na.rm=T)==TRUE])
  variablesinitialNAs <- length(apply(apply(datos,2,function(x) is.na(x)),2,any)[apply(apply(datos,2,function(x) is.na(x)),2,any)==TRUE])
  datos[datos==0] <- NA
  matriz2 <-matriz
  matriz2[initialrow:nrow(matriz),initialcolumn:ncol(matriz)] <-datos 
  rm(datos) 
  #### Cero - NA Split ####
  if(treatment==1) vectortratamientos <- (matriz[initialrow:nrow(matriz),treatment1col])
  if(treatment==2) vectortratamientos <- (matriz[initialrow:nrow(matriz),treatment2col])
  if(treatment==3) vectortratamientos <- (paste(matriz[initialrow:nrow(matriz),treatment1col],matriz[initialrow:nrow(matriz),treatment2col], sep="&"))
  
  vectortratamientos<-factor(vectortratamientos,levels = unique(vectortratamientos))
  
  datos <- split(as.data.frame(matriz2),vectortratamientos,drop=T) 
  processedlist <- lapply(datos, function(x) NAorZeroTreatment(x,threshold,initialcolumn))
  matriz2 <- do.call("rbind",  unname(processedlist))
  matriz2 <- matriz2[match(rownames(matriz), rownames(matriz2)), ]
  
  finalNAs <-sum(is.na(matriz2[initialrow:nrow(matriz2),1:ncol(matriz2)]))
  final0s <- sum(matriz2[initialrow:nrow(matriz2),1:ncol(matriz2)]==0,na.rm=TRUE)
  variablesfinalNAs <- length(apply(apply(matriz2,2,function(x) is.na(x)),2,any)[apply(apply(matriz2,2,function(x) is.na(x)),2,any)==TRUE])
  if(finalNAs!=0){
    cat("\nOriginal data contained", initial0s, "zeroes and", initialNAs,"NAs in ",variablesinitial0s ," and ",variablesinitialNAs," variables, respectively. After processing", finalNAs, "values present in ", variablesfinalNAs, "variables have been considered suitable for imputation according to the defined",threshold,"threshold.\n")
  }
  else{
    cat("\nOriginal data contained", initial0s, "zeroes and", initialNAs,"NAs in ",variablesinitial0s ," and ",variablesinitialNAs," variables, respectively. After processing no variable has been considered suitable for imputation according to the defined",threshold,"threshold.\n")  
    return(matriz2)
  }
  
  #### Imputacion ####
  
  if(imputation=="none"){
    matriz2[is.na(matriz2)]<-0
    return(matriz2)
  } 
  else if(imputation=="KNN"){
    #suppressMessages(suppressWarnings(require(impute)))
    imputed<-matriz2[initialrow:nrow(matriz2),initialcolumn:ncol(matriz2)]
    imputed<-as.data.frame(impute::impute.knn(as.matrix(imputed),k)$data)
    
    
  } 
  else if(imputation=="RF") {
    #suppressMessages(suppressWarnings(require(missForest)))
    imputed<-matriz2[initialrow:nrow(matriz2),initialcolumn:ncol(matriz2)]
    if(cores>=2){
      #suppressMessages(suppressWarnings(require(doParallel)))
      doParallel::registerDoParallel(cores=cores) 
      #suppressMessages(suppressWarnings(require(misForest)))
      imputed<-(missForestedit(imputed,parallelize="variable"))$ximp
      #matriz2[initialrow:nrow(matriz2),initialcolumn:ncol(matriz2)]<-imputed
    } else if(cores==1){
      imputed<-(missForestedit(imputed,parallelize="no"))$ximp}
    #matriz2[initialrow:nrow(matriz2),initialcolumn:ncol(matriz2)]<-imputed
  } 
  else{
    return(cat("\nError. Please revise selected imputation method."))
  }
  if(any(imputed<0)&all(matriz[initialrow:nrow(matriz),initialcolumn:ncol(matriz)]>=0)){
    imputed[imputed<0]<-0
  }
  matriz2[initialrow:nrow(matriz2),initialcolumn:ncol(matriz2)]<-imputed
  cat("Missing values have been imputed according to",imputation,"method.\n")
  return(matriz2)
}
