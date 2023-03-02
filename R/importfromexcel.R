#
# Copyright 2020 Luis Valledor / Laura Lamelas
# Last revision 25.11.2020
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

#' @name importfromexcel
#' @title Import datasets to R from an Excel file
#' @description This function is used to import and check datasets from Excel files
#' @usage importfromexcel("datasetname")
#' @param datasetname Character. The name assigned to the imported data
#' @details This function imports a spreadsheet from office or open office (xlsx or xls extensions) to R,
#' creates a dataset list to be further used in pRocessomics.
#' It also performs some initial quality and consistency checks. If data is not OK
#' it will give the user some hints for correcting it.
#' The different omic levels within your dataset (from 1 to n) must be integrated
#' in the same Excel file, one sheet per level. Each Excel sheet should be named
#' with the name you want to give to each omic level.
#' As a golden rule: cases (samples) data should be in rows, while variables should
#' be in columns. First column must indicate the name of the cases or samples, which
#' will be further treated as rownames in the rest of the functions of this package.
#' Sample names should be the same across the different omic levels. All omic levels
#' should have the same number of cases (this is a current limitation of this package).
#' We recommend that second and third column will correspond to the different treatments
#' you applied to your data, these factors must be also maintained all across the
#' different omic levels. After these columns, variables (proteins,
#' metabolites, rnas, ...) should be indicated.
#' We strongly recommend that you employ unique protein, metabolite, gene, ... accessions
#' instead of definitions or descriptions for naming variables. Later you can upload
#' an annotation file to get the most of your data with \code{\link{importannotation}} function.
#' Please go to \url{http://processomics.valledor.info} or \url{https://github.com/Valledor/pRocessomics/}
#' to download a template you can use to fill with your data.
#' @return pRoDS class list containing all omic levels of the dataset and a vector with the dataset structure.
#' @author Luis Valledor and Laura Lamelas
#' @export
#' @importFrom utils read.delim select.list 
#' @importFrom readxl excel_sheets read_excel
#' @seealso importannotation




importfromexcel <- function(datasetname = NULL) {
  #### Starting console logging ####
  diahora <-gsub(":","",gsub(" ","_", gsub("-","",Sys.time())))
  diahora <-substr(diahora,3,nchar(diahora)-2)
  sinkfilename <- paste(diahora,"_ImportLog.txt",sep="")
  sink(sinkfilename, split = T)
  
  headings_wizards("Welcome to pRocessomics Data Import Wizard")
  texts_wizard("\n\nThis wizard is aimed to guide you for properly import your data from excel files.")
  texts_wizard(
    "\nPlease note that excel data should start at A1 cell in every book sheet. In addition, decimal separators must a point (.), not comma (,)."
  )
  texts_wizard("\nPress [enter] to continue...")
  invisible(readline())
  if (is.null(datasetname) == TRUE) {
    texts_wizard("Please provide a name for your dataset (not numeric: 123 -> wrong; 123a -> right) and press [enter]:")
    datasetname <- readline("")
  }
  ### no need to quote the name, readline() forces as.character
  if (datasetname == "")
    stop(
      "Please provide a name for your dataset (i.e. mydataset)")
  if (is.character(datasetname) != TRUE)
    stop(
      "Please provide a character-class name for your dataset (i.e. mydataset)")
  if (is.na(datasetname) == TRUE)
    stop(
      "Please provide a name for your dataset (i.e. mydataset).")
  if (is.numeric(datasetname) == TRUE)
    stop("Don't use numbers, please provide a character name for your dataset (i.e. mydataset).")
  
  if (exists(datasetname, envir = .GlobalEnv) == TRUE) {
    texto <-
      paste(
        "\n",
        datasetname,
        " object already exists in global environment and it will be over-written with you new data imput. Do you want to continue?"
      )
    texts_wizard(texto)
    answer <- utils::select.list(choices = truefalse)
    if (answer == "Yes")
      answer <- TRUE
    if (answer == "No")
      answer <- FALSE
    continue <- answer
    if (continue == F)
      stop("Import aborted by user")
  }
  
  texts_wizard("\nNow select the xls or xlsx file containg your dataset(s).\n\n")
  readline(prompt = "Press [enter] to open the window (note that it can be opened behind this window).")
  excelfilename <- file.choose()
  excelfilename <- file.path(excelfilename)
  
  
  fileextension <- strsplit(excelfilename, "[.]")[[1]]
  fileextension <- fileextension[[length(fileextension)]]
  fileextension <- toupper(fileextension)
  if (fileextension %in% c("XLS", "XLSX") == F) {
    stop("\n\nIMPORT ERROR: File to import must be in xlsx or xls formats.\n\n")
  }
  if (!file.exists(excelfilename)) {
    stop(
      "\n\nIMPORT ERROR: File doesn't exist. Please provide a valid filename to import.\nHint: Special characters in path, folder or filename may be a problem... "
    )
  }
  #### Imports from excel ####
  sheetnames <- readxl::excel_sheets(excelfilename) #get sheet names
  dataset <-
    lapply(sheetnames, function(x)
      readxl::read_excel(excelfilename, sheet = x, .name_repair = "unique"))
  dataset <- lapply(dataset, as.data.frame)
  
  
  #### Select columns, rows, and treatments ####
  # We ask the user to input starting row, starting column, treatment 1 and treatment 2.
  # We provide a dataset head in order to help user to choose correct columns.
  texts_wizard(
    "\nBelow you will find a small portion of your dataset so you can easily check which is your starting row, starting column, and treatment columns.\n"
  )
  rows_name <- c(" ", "Row 1", "Row 2", "Row 3", "Row 4", "Row 5")
  cols_name <- c("Col 1", "Col 2", "Col 3", "Col 4", "Col 5")
  if(length(which(dim(dataset[[1]]) <= 5)) == 0){
    datasetsample <- dataset[[1]][1:5, 1:5]
  }
  if(length(which(dim(dataset[[1]]) < 5)) %in% c(1,2)){
    if(nrow(dataset[[1]]) > 5) {
      datasetsample <- dataset[[1]][1:5,1:ncol(dataset[[1]])]
    }
    else if(ncol(dataset[[1]]) > 5) {
      datasetsample<- dataset[[1]][1:nrow(dataset[[1]]),1:5]}
    else{
      datasetsample <- dataset[[1]][1:nrow(dataset[[1]]), 1:ncol(dataset[[1]])]
    }
  }
    datasetsample <- rbind(colnames(datasetsample), datasetsample)
    rownames(datasetsample) <- rows_name[1:nrow(datasetsample)]
    colnames(datasetsample) <- cols_name[1:ncol(datasetsample)]
  
  print(datasetsample)
  texts_wizard("\nEnter the number of the column containing sample names and press [enter]: ")
  samplenamescolumn <- as.integer(readline())
  checkuserinput(samplenamescolumn, nrow = NULL, ncol(dataset[[1]]))
  texts_wizard("\nEnter the number of the column with the first treatment and press [enter]: ")
  treatment1col <- as.integer(readline())
  checkuserinput(treatment1col , nrow = NULL, ncol(dataset[[1]]))
  texts_wizard(
    "\nEnter the number of the column with the second treatment and press [enter] (optional, press [enter] if you only have one treatment/level): "
  )
  treatment2col <-  readline()
  if (treatment2col != "") {
    treatment2col <- as.integer(treatment2col)
    checkuserinput <-
      checkuserinput(treatment2col, nrow = NULL, ncol(dataset[[1]]))
    if (treatment2col == treatment1col)
      treatment2col = NULL
  }
  if (treatment2col == "")
    treatment2col <- NULL
  texts_wizard(
    "\nEnter the number of the column with the first variable (numerical, not treatment) and press [enter]:  "
  )
  initialcolumn <- as.integer(readline())
  checkuserinput(initialcolumn , nrow = NULL, ncol(dataset[[1]]))
  texts_wizard("\nEnter the number of the row with the first case (individual) and press [enter]: ")
  initialrow <- as.integer(readline())
  checkuserinput(initialrow, nrow(dataset[[1]]), ncol = NULL)
  
  
  if (is.null(treatment1col))
    stop("\n\nIMPORT ERROR: Please select at least one column indicating treatments.\n\n")
  
  if (is.null(initialcolumn)) {
    
    texts_wizard("\nTABLE ERROR: You should indicate the first numeric data column, argument initialcol.\n")

    texts_wizard(
      "\nHint: To ease things, you should provide data in a way pRocessomics can understand. Please check table organization in the documentation section. Data block should start in one column and take all columns to the end of the table. Initial columns must define Individual/sample name (column 1), treatments (columns 2 and 3), some atributes (n cols) and then start with data.\n\n"
    )
    stop()
  }
  
  #### Remove all rows but initial ####
  if (is.null(initialrow)) {
    texts_wizard(
      "\n\nWARNING: initial row was not provided. Considering initialrow=1 (Data of first case is in row 1). Please cancel and set proper initialrow if you don't like this setting.\n\n"
    )
    readline(prompt = "Press [enter] to continue")
    initialrow = 1
  }
  if (initialrow != 1) {
    texts_wizard(
      "\n\nRemoving rows without any relevant data for downstream analyses (rows 1:initialrow).\n"
    )
    dataset <-
      lapply(dataset, function(x) {
        x <- x[initialrow:nrow(x), ]
        return(x)
      })
    initialrow = 1
  }
  
  #### Check if sample names are equal in all datasets ####
  samplenamesequal <- arecolumnsequal(dataset, 1)
  if (all(samplenamesequal) == FALSE) {
    a <- as.character(paste(
      "\nTABLE ERROR: Error in case names in",
      sheetnames[samplenamesequal == F],
      "dataset(s).\n"
    ))
    texts_wizard(a)
    texts_wizard(
      "\nHint: Case (sample) names should be indicated in the first column of each tables, futhermore these names should be the same across datasets. Cases should also have the same order in all sheets.\n\n"
    )
    stop()
  }
  
  #### Check if treatment 1 names are equal in all datasets ####
  if (length(dataset) > 1) {
    #Only do that if the number of elements of the list is >=2
    treat1equal <- arecolumnsequal(dataset, treatment1col)
    if (all(treat1equal) == FALSE) {
      a <-
        as.character(paste(
          "\nTABLE ERROR: Treatment 1 column is not consistent among different datasets. Check names in ",
          sheetnames[treat1equal == F],
          "dataset(s).\n"
        ))
      texts_wizard(a)
      texts_wizard(
        "\nHint: Treatment names should be the same across datasets. Samples should also have the same order in all sheets.\n Please correct erros and re-run.\n\n"
      )
      stop()
    }
    #### Check if treatment 2 names are equal in all datasets ####
    if (!is.null(treatment2col)) {
      treat2equal <- arecolumnsequal(dataset, treatment2col)
      if (all(treat2equal) == FALSE) {
        a <-
          as.character(paste(
            "\nTABLE ERROR: Treatment 2 column is not consistent among different datasets. Check names in ",
            sheetnames[treat2equal == F],
            "dataset(s).\n"
          ))
        texts_wizard(a)
        texts_wizard(
          "\nHint: Treatment names should be the same across datasets. Samples should also have the same order in all sheets.\n Please correct errors and re-run.\n\n"
        )
        stop()
      }
    }
  }
  
  
  #
  if (length(dataset[[1]][, samplenamescolumn]) != length(unique(dataset[[1]][, samplenamescolumn]))) {
    texts_wizard(
      "\nTABLE ERROR: Sample names must be unique. Sample names will be renamed\n"
    )
    samplenames <-
      make.names(dataset[[1]][, samplenamescolumn], unique = T)
    dataset <-
      lapply(dataset, function(x) {
        rownames(x) <- as.vector(samplenames) #sample names as rownames
        return(x)
      })
    
  } else {
    dataset <-
      lapply(dataset, function(x) {
        rownames(x) <-
          as.vector(x[, samplenamescolumn]) #sample names as rownames
        return(x)
      })
  }
  
  dataset <-
    lapply(dataset, function(x) {
      x <- x[, -1] #removal of samplenames column
      return(x)
    })
  
  
  #### Removing Empty columns ####
  texts_wizard("\n\nRemoving empty columns...\n")
  names(dataset) <- sheetnames
  dataset <-
    lapply(seq_along(dataset), function(x)
      RemoveEmptyColumns_import(dataset[[x]], initialrow, initialcolumn, names(dataset)[[x]]))
  
  #### Check empty rows ####
  
  emptyrows <-
    lapply(dataset, function(x)
      CheckEmptyRows(x, initialrow, initialcolumn) == TRUE)
  emptyrows <- unlist(emptyrows)
  if (length(emptyrows[emptyrows == T]) > 0) {
    
      a <- as.character(paste(
      "\nFATAL ERROR: There are some cases without any numeric data: ",
      paste(names(emptyrows[emptyrows == T]), collapse = ", "),
      "\n"
    ))
      texts_wizard(a)
    texts_wizard(
      "\nHint: Please make sure that your dataset is complete. In at least one case and treatment, all of your variables don't have any significant numerical value.\n\n"
    )
    stop()
  }
  
  #### Is data numeric? We must consider numbers and NA ####
  mydatanumeric <-
    lapply(dataset, function(x)
      isnumeric(x, initialcolumn))
  mydatanumeric <- unlist(mydatanumeric)
  if (all(mydatanumeric) == FALSE) {
    
    a <- 
      as.character(paste(
        "\nTABLE ERROR: There are some numeric variables that are not numeric. Please check the following variables across your dataset:\n\n"
      ))
    texts_wizard(a)
    cat(names(mydatanumeric)[mydatanumeric == F])
    texts_wizard(
      "\n\nHint: Make sure that the decimal separator is set properly. Be careful since R employs period (.) as decimal separator.\n\n"
    )
    stop()
  }
  names(dataset) <- sheetnames
  if (is.null(treatment2col) == TRUE)
    treatment2col <- NA
  defaultsettings <-
    c(
      initialrow = initialrow,
      initialcolumn = initialcolumn - 1,
      treatment1col = treatment1col - 1,
      treatment2col = treatment2col - 1,
      treatment = NA
    )
  class(dataset) <- "pRoDS"
  class(defaultsettings) <- "dsettings"
  
  exportlist <- list(dataset, defaultsettings)
  names(exportlist) <-
    c(datasetname, paste(datasetname, "_dsettings", sep = ""))
  list2env(exportlist, envir = .GlobalEnv)
  
  #### Exporting session log ####
  texts_wizard("\nimportfromexcel: Duty Accomplished.")
  dog <- as.character(paste("Session log, ",sinkfilename, " has been saved in your working directory: ", getwd(),sep=""))
  texts_wizard(dog)
  sink()
  
}
