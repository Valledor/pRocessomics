#' @name export_table
#' @title Function for table exportation
#' @description This function exports pRocessomics generated tables to Excel format
#' @usage export_table(result_list,filename=NULL)
#' @param result_list name of the list to be exported in the local environment
#' @param filename quoted name to be given to the result_list object with the extension .xlsx
#' @author Luis valledor and Laura Lamelas
#' @export

#' @importFrom writexl write_xlsx
#' @importFrom mixOmics plotVar
#' @importFrom grDevices pdf dev.off 

export_table <- function(result_list,filename=NULL){
  #suppressPackageStartupMessages(library(writexl))
  if(inherits(result_list, "POL")){
    texts_wizard("\nEXPORTING PRE-PROCESSED, TRANSFORMED and/or FEATURE SELECTED TABLES")
    texts_wizard("\nA single core will be used. It may take a while...\n ")
    for(i in 1:length(result_list)){
      texts_wizard(paste("Generating and saving", names(result_list)[i], "table. Data was grouped.\n"))
      if(is.null(filename)) filename= "preprocessed_table.xlsx"
      result_list2<-lapply(result_list, function(x) as.data.frame(x))
      result_list2<-lapply(result_list2, function(x) as.data.frame(cbind("Samples"=rownames(x),x)))
      writexl::write_xlsx(result_list2, path=filename)
      texts_wizard("Done!\n")
    }
  }
  
  if(inherits(result_list, "vennanalysis")){
    texts_wizard("\nEXPORTING VENN TABLES")
    texts_wizard("\nA single core will be used. It may take a while...\n ")
    texts_wizard("\nGenerating and saving Venn diagram associated table.\n")
    if(is.null(filename)) filename= "Venn_table.xlsx"
    result_list2<-as.data.frame((result_list$data))
    writexl::write_xlsx(result_list2, path=filename)
    cat("Done!\n")
  }
  
  if(inherits(result_list, "UNIVAR")){
    result_list <- result_list$composite
    result_list2<-lapply(result_list, function(x) as.data.frame(x))
    result_list2<-lapply(result_list2, function(x) as.data.frame(cbind("Variables"=rownames(x),x)))
    
    texts_wizard("EXPORTING UNIVARIABLE TABLES")
    texts_wizard("\nA single core will be used. It may take a while...\n")
    
    for(i in 1:length(result_list)){
      texts_wizard(paste("Generating and saving", names(result_list)[i], "table. Data was grouped. Mean, SD, and p and q-values will be provided.\n"))
      if(is.null(filename)) filename= "univar_table.xlsx"
      writexl::write_xlsx(result_list2, path=filename)
      
      texts_wizard("Done!\n")
      }
    }
      
  if(inherits(result_list, "MMOG")){
    cat("EXPORTING MAPMAN TABLES\nA single core will be used. It may take a while...\n")
    cat("Generating and saving Mapman table. Data was grouped. Mean and SD will be provided.\n")
    
    colnames(result_list$mean) <-paste(colnames(result_list$mean),"-Mean",sep="")
    colnames(result_list$sd) <-paste(colnames(result_list$sd),"-SD",sep="")
    if(length(result_list)==5){
      pvals<-as.data.frame(result_list$pvalues)
      colnames(pvals) <-c("p-values")
    }
    
    newtable <-cbind(result_list$mean[1],result_list$sd[1])
    if(ncol(result_list$mean)>1){ #esto va a ser length me da a mi
      for(i in 2:ncol(result_list$mean)){
        newtable <-cbind(newtable,result_list$mean[i],result_list$sd[i])
      }
    }
    if(is.null(filename)) filename= "mapman_table.xlsx"
    #write.xlsx(newtable, file=filename, sheetName="MapmanGroups")
    cat("Done!\n")
    }
  
  if(inherits(result_list, "MMO")){
    cat("EXPORTING MAPMAN TABLES\nA single core will be used. It may take a while...\n ")
    cat("Generating and saving Mapman table. Data was ungrouped by treatment. Mean will be provided.\n")
    colnames(result_list$mean) <-paste(colnames(result_list$mean),"-Mean",sep="")
    if(is.null(filename)) filename= "mapman_table.xlsx"
    Mapman<-cbind(rownames(result_list$mean),result_list$mean)
    writexl::write_xlsx(Mapman, path=filename)
    cat("Done!\n")
  }
  
  if(inherits(result_list, "pcaanalysis")){
    cat("EXPORTING PCA TABLES\nA single core will be used. It may take a while...\n ")
    if(is.null(filename)) filename= "pca_table.xlsx"
    exportpcatable<-rbind(round((result_list$pca$sdev^2)*100/sum(result_list$pca$sdev^2),2),cumsum(round((result_list$pca$sdev^2)*100/sum(result_list$pca$sdev^2),2)),result_list$pca$sdev)
    row.names(exportpcatable)<-c("Proportion_of_Variance","Cumulative_Proportion_of_Variance","Standard_Deviation")
    colnames(exportpcatable)<-paste("PC",c(1:length(result_list$pca$sdev)),sep="")
    result_list <- list("Scores"=result_list$pca$x, "Loadings"=result_list$pca$rotation, "Explained Variance"=exportpcatable)
    result_list2<-lapply(result_list, function(x) as.data.frame(x))
    result_list2<-lapply(result_list2, function(x) as.data.frame(cbind(" "=rownames(x),x)))
    writexl::write_xlsx(result_list2, path=filename)
    cat("Generating and saving Scores table. Step 1/3\n")  
    cat("Generating and saving Loadings table. Step 2/3\n")
    cat("Generating and saving Explained Variance table. Step 3/3\n")
    cat("Done!\n")
  }
  
  if(inherits(result_list, "icaanalysis")){
    cat("EXPORTING ICA TABLES\nA single core will be used. It may take a while...\n ")
    if(is.null(filename)) filename= "ica_table.xlsx"
    exporticatable<-data.frame(ExplVar=round(result_list$ica$vafs*100,2), Cumulat=round(cumsum(result_list$ica$vafs*100),2)) 
    rownames(exporticatable)<-paste0("IC",1:nrow(exporticatable))
    result_list2 <- list("Signal Estimates"=result_list$ica$S, "Unmixing matrix"= t(result_list$ica$W),"Variance accounted"=exporticatable )  
    cat("Generating and saving signal estimates table. Step 1/3\n")
    cat("Generating and saving estimated unmixing matrix. Step 2/3\n")
    cat("Generating and saving variance-accounted-for by each component. Step 3/3\n")  
    
    result_list2<-lapply(result_list2, function(x) as.data.frame(x))
    result_list2<-lapply(result_list2, function(x) as.data.frame(cbind(" "=rownames(x),x)))
    writexl::write_xlsx(result_list2, path=filename)
    cat("Done!\n")
    
  }
  
  if(inherits(result_list, "mcoaanalysis")){
    cat("EXPORTING MCOA TABLES\nA single core will be used. It may take a while...\n ")
    if(is.null(filename)) filename= "mcoa_table.xlsx"
    norm_eigenvalues <- result_list$mcoa$mcoa$pseudoeig/sum(result_list$mcoa$mcoa$pseudoeig) #Sacamos los autovalores, y os dividimos por la suma para tener el porcentaje de varianza en teoria recogida por cada tratamiento
    screeplotdata <- data.frame(Cs=paste("Comp",c(1:length(norm_eigenvalues)),sep=" "), ExplVar=norm_eigenvalues, Cumulat=cumsum(norm_eigenvalues)*100/sum(norm_eigenvalues),stringsAsFactors = FALSE) #Creamos tabla de datos con % varianza y % var acumulado
    screeplotdata$Cs <- factor(screeplotdata$Cs, levels = unique(screeplotdata$Cs)[order(screeplotdata$ExplVar, decreasing = TRUE)]) # Ordenamos los PCs segun valor
    rownames(screeplotdata)<-screeplotdata[,1]
    screeplotdata<-screeplotdata[,-1]
    result_list2<-list("SynVar Scores"=result_list$mcoa$mcoa$SynVar,"Co-inertia Norm Score"=result_list$mcoa$mcoa$Tl1,
                       "Coordinates onto Syn Scores"=result_list$mcoa$mcoa$Tco,"Pseudo-eigenvalues"=result_list$mcoa$mcoa$cov2,
                       "Cumulative pseudo-eigenvalues"=screeplotdata)
    
    result_list2<-lapply(result_list2, function(x) as.data.frame(x))
    result_list2<-lapply(result_list2, function(x) as.data.frame(cbind(" "=rownames(x),x)))
    writexl::write_xlsx(result_list2, path=filename)
    
    cat("Generating and saving synthetic scores table. Step 1/5\n")
    cat("Generating and saving co-inertia normed scores. Step 2/5\n")
    cat("Generating and saving column coordinates onto synthetic scores. Step 3/5\n")
    cat("Generating and saving pseudo-eigenvalues. Step 4/5\n")
    cat("Generating and saving cumulative pseudo-eigenvalues. Step 5/5\n")
    cat("Done!\n")
  }
  
  if(inherits(result_list, "kmeansanalysis")){
    cat("EXPORTING K-MEANS GROUPING TABLE\nA single core will be used. It may take a while...\n ")
    if(is.null(filename)) filename= "kmeans_table.xlsx"
    cat("Saving xlsx file. Step 1/1\n")
    writexl::write_xlsx(x = result_list$kmeans_list,path = filename)
  }
  
  if(inherits(result_list, "splsanalysis")){
    spls_object<-result_list[[1]]
    cat("EXPORTING SPLS TABLES\nA single core will be used. It may take a while...\n ")
    cat("Generating and saving scores tables. Step 1/3\n")
    if(is.null(filename)) filename= "spls_table.xlsx"
    grDevices::pdf(file = NULL)
    vartable <-invisible(mixOmics::plotVar(spls_object))
    grDevices::dev.off()
    vartable<- vartable[,1:3]
    Expl<-as.data.frame(unlist(spls_object$explained_variance))
    colnames(Expl)<-c("ExplainedVariance")
    result_list2<-list("Scores-X"=spls_object$variates$X,"Scores-Y"=spls_object$variates$Y,
                       "Loadings-X"=spls_object$loadings$X,"Loadings-Y"=spls_object$loadings$Y,
                       "VarScores"=vartable,"Expl.Var"=Expl)
    result_list2<-lapply(result_list2, function (x) as.data.frame(x))
    result_list2<-lapply(result_list2, function (x) as.data.frame(cbind(" "=rownames(x),x)))
    cat("Generating and saving variable loadings and scores tables. Step 2/3\n")
    cat("Generating and saving Explained Variance table. Step 3/3\n")
    writexl::write_xlsx(result_list2, path=filename)
  }
  
  if(inherits(result_list, "daanalysis")){
    da_object<-result_list[[1]]
    cat("EXPORTING DA TABLES\nA single core will be used. It may take a while...\n ")
    cat("Generating and saving scores tables. Step 1/3\n")
    if(is.null(filename)) filename= "da_table.xlsx"
    grDevices::pdf(file = NULL)
    vartable <-invisible(mixOmics::plotVar(da_object))
    grDevices::dev.off()
    vartable<- vartable[,1:3]
    Expl<-as.data.frame(unlist(da_object$prop_expl_var))
    colnames(Expl)<-c("ExplainedVariance")
    result_list2<-list("Scores-X"=da_object$variates$X,"Scores-Y"=da_object$variates$Y,
                       "Loadings-X"=da_object$loadings$X,"Loadings-Y"=da_object$loadings$Y,
                       "VarScores"=vartable,"Expl.Var"=Expl)
    result_list2<-lapply(result_list2, function (x) as.data.frame(x))
    result_list2<-lapply(result_list2, function (x) as.data.frame(cbind(" "=rownames(x),x)))
    cat("Generating and saving variable loadings and scores tables. Step 2/3\n")
    cat("Generating and saving Explained Variance table. Step 3/3\n")
    writexl::write_xlsx(result_list2, path=filename)
  }
  
  if(inherits(result_list, "bdaanalysis")){
    diablo_object<-result_list[[1]]
    cat("EXPORTING DIABLO ANALYSIS TABLES\nA single core will be used. It may take a while...\n ")
    cat("Generating and saving scores tables. Step 1/3\n")
    if(is.null(filename)) filename= "diablo_table.xlsx"
    vars<-diablo_object$variates
    names(vars)<-paste0("Score-",names(vars),sep="")
    cat("Generating and saving variable loadings and scores tables. Step 2/3\n")
    loads<-diablo_object$loadings
    names(loads)<-paste0("Loadings-",names(loads))
    grDevices::pdf(file = NULL)
    vartable <-(mixOmics::plotVar(diablo_object))
    grDevices::dev.off()
    vartable<-list("VarScores"=vartable[,1:3])
    cat("Generating and saving Explained Variance table. Step 3/3\n")
    expvar<-as.data.frame(unlist(diablo_object$explained_variance))
    colnames(expvar)<-c("Explained.Variance")
    result_list2<-c(vars,loads,vartable,expvar)
    result_list2<-lapply(result_list2, function (x) as.data.frame(x))
    result_list2<-lapply(result_list2, function (x) as.data.frame(cbind(" "=rownames(x),x)))
    writexl::write_xlsx(result_list2,path=filename)
  }
  
  if(inherits(result_list, "parametricity")){
    if(is.null(filename)) filename= "parametricity_table.xlsx"
    result_list2 <- mapply(cbind, result_list$normal.distribution, result_list$homocedasticity, SIMPLIFY = F)
    result_list2<-lapply(result_list2, function(x) as.data.frame(cbind("Variable.IDs"=rownames(x),x)))
    writexl::write_xlsx(result_list2, path=filename)
  }
}

