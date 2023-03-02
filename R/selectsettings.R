#
# Copyright 2019 Luis Valledor / Laura Lamelas
# Last revision 31.03.2020
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

selectsettings <- function(originaldatasetname,datalist){
  options(warn = -1)
  originaldatasetname<-strsplit(originaldatasetname,"_")[[1]][1]
#### Setting file exists
  if(exists(paste(originaldatasetname,"_dsettings",sep=""), envir = .GlobalEnv)==TRUE){
    settings <- get(paste(originaldatasetname,"_dsettings",sep=""),envir = .GlobalEnv)
    texts_wizard("\nImporting dataset settings from global environment")
    dog<-as.character(paste("\nInitialrow:", settings[1], "\tTreatment1 column: ", settings[3], "\tTreatment2:",settings[4],"\tInitialcolumn: ",  settings[2], sep=" "))
    texts_wizard(dog)
    dog<-as.character(paste("\nSelected Treatment:",settings[5], sep=" "))
    texts_wizard(dog)
    texts_wizard("\nDo you agree employing these dataset settings?")
    truefalse <- c("Yes","No")
    answ <- utils::select.list(choices = truefalse)
    if(answ == "Yes") answ <- TRUE
    if(answ == "No") answ <- FALSE
    
    if(answ==T){
      return(settings)
    } 
    else {
      if(is.na(settings[4])==TRUE|settings[3] == settings[4]){
        treatment <- 1
        settings[5] <- 1
        texts_wizard("\nOnly one treatment column was defined. Exporting previous settings. Press [enter] to continue.")
        exportname <-paste(originaldatasetname,"_dsettings",sep="")
        exportlist <- list(settings)
        names(exportlist) <- exportname
        list2env(exportlist,envir=.GlobalEnv)
        invisible(readline())
        return(settings)
      }
      if(!is.na(settings[4])){
        treatmentlist <- c("Treatment 1", "Treatment 2", "Combination of both")
        treatment <-  utils::select.list("Please select treatment to be used: ",choices = treatmentlist,multiple = F,graphics = F)  
        treatment <- which(treatmentlist==treatment)
        settings[5] <- treatment
        texts_wizard("\nDo you want to save these settings for other analyses?")
        answer <- utils::select.list(choices = truefalse)
        if(answer == "Yes"){
          exportname <-paste(originaldatasetname,"_dsettings",sep="")
          exportlist <- list(settings)
          names(exportlist) <- exportname
          list2env(exportlist,envir=.GlobalEnv)
          }
        return(settings)
      }
      
    } 
    
    return(settings)
  }

  #### Setting file does not exist
  
  texts_wizard("\n\nInformation about dataset structure is not available.")
  texts_wizard("\nDo you want to introduce these settings by hand?")
  answer <- utils::select.list(choices = truefalse)
  if(answer == "No" ) stop("Analysis aborted by user")
  
  texts_wizard("\nBelow you will find a small portion of your dataset so you can easily check which is your starting row, starting column, and treatment columns.\n")
  datasetsample<-datalist[[1]][1:5,1:5]
  datasetsample<-rbind(colnames(datasetsample),datasetsample)
  datasetsample<-cbind(rownames(datasetsample),datasetsample)
  rownames(datasetsample) <- c(" ", "Row 1","Row 2","Row 3","Row 4","Row 5")
  colnames(datasetsample) <- c(" ", "Col 1","Col 2","Col 3","Col 4","Col 5")
  print(datasetsample)
  texts_wizard("Enter the number of the column with the first treatment and press [enter]: ")
  treatment1col <- as.integer(readline())
  checkuserinput(treatment1col ,nrow=NULL,ncol(datalist[[1]]))
  texts_wizard ("Enter the number of the column with the second treatment and press [enter]: (optional, press [enter] if you only have one treatment/level): ")
  treatment2col <-  readline()
  
  if(treatment2col!="") {
    treatment2col <- as.integer(treatment2col)
    checkuserinput <-
      checkuserinput(treatment2col, nrow = NULL, ncol(datalist[[1]]))
    if (treatment2col == treatment1col) {
      treatment2col = NULL
    }
  }
 
  if(treatment2col=="")  treatment2col <- NULL
  if(is.null(treatment2col)) {
    treatment=1
    treatment2col=NA
    }
  
  texts_wizard("Enter the number of the column with the first variable (numerical, not treatment) and press [enter]:  ")  
  initialcolumn <- as.integer(readline())
  checkuserinput(initialcolumn ,nrow=NULL,ncol(datalist[[1]]))
  texts_wizard("Enter the number of the row with the first case (individual) and press [enter]: ")
  initialrow <- as.integer(readline())
  checkuserinput(initialrow,nrow(datalist[[1]]),ncol=NULL)
  #print(c(initialrow,initialcolumn,treatment1col,treatment2col))
  if(!is.na(treatment2col)){
    treatmentlist <- c("Treatment 1", "Treatment 2", "Combination of both")
    texts_wizard("Please select treatment to be used: ")
    treatment <-  utils::select.list(choices = treatmentlist,multiple = F,graphics = F)  
    treatment <- which(treatmentlist==treatment)
    settings <- c(initialrow, initialcolumn, treatment1col, treatment2col, treatment)
    texts_wizard("\nDo you want to save these settings for other analyses?")
    answer <- utils::select.list(choices = truefalse)
    if(answer == "Yes"){
      exportname <-paste(originaldatasetname,"_dsettings",sep="")
      exportlist <- list(settings)
      names(exportlist) <- exportname
      list2env(exportlist,envir=.GlobalEnv)
    }
    return(settings)
  }
  settings <- c(initialrow, initialcolumn, treatment1col, treatment2col, treatment)
  texts_wizard("\nDo you want to save these settings for other analyses?")
  answer <- utils::select.list(choices = truefalse)
  if(answer == "Yes"){
    exportname <-paste(originaldatasetname,"_dsettings",sep="")
    exportlist <- list(settings)
    names(exportlist) <- exportname
    list2env(exportlist,envir=.GlobalEnv)
  }
  return(settings)
  

}
