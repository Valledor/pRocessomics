#### To DOs ####

# Cambiar cores=2 por parallel? En ese caso poner max number of cores
# De no ser asi, y cores es dado por usuario, si cores >numero cores, reducirlo a este numero
    #esto te lo hace el paquete solo, bueno por lo que he leido te da el mensaje de error, coge los cores con:
    #getDoParWorkers()

# IQR, hacerlo con medias de tratamiento es un poco absurdo, no? cambiarlo para tomar todos los valores?
# IQRfilter, treatment=0 hace esto

# Las individuales funcionan, falta testar preprocess_omic_list, transfsel_omic_list, featsel_omic_list

# Seleccionar dataset que se use de ejemplo a lo largo del paquete
    #estaria bien tener uno de prueba con varias capas omicas... 

# Documentar, documentar, documentar

# Tests unitarios, complementar? dar por buenos?


preprocess_omic_list <- function(datalist, initialrow=1, initialcolumn=2, treatment1col=1, treatment2col=2, treatment=1, imputation=NULL, abdbal=NULL, threshold=0.25, k=3, parallel=FALSE){
  #### Initial checks ####
  # Dataset names, number of rows, cases, etc.
  if(hasArg(datalist)==FALSE) stop("\nPlease introduce a valid dataset\n")
  if(class(datalist) != "list")
    stop("A list of matrices or dataframes corresponding to each level is expected")
  if(is.null(names(datalist)))
    stop("Names of elements in list are missing. Please provide names for each level")
  if(all(sapply(lapply(datalist, function(x) t(x[initialrow:nrow(x),initialcolumn:ncol(x)])), dim)[2,]==sapply(lapply(datalist, function(x) t(x[initialrow:nrow(x),initialcolumn:ncol(x)])), dim)[2,1])==FALSE)
    stop("The different matrices have an unequal number of individuals. The number of cases of each level should be the same. Please read instructions for more information")
  dataframesnames<-lapply(datalist, function(x) rownames(x))
  if(all(unlist(lapply(dataframesnames, function(x) identical(dataframesnames[[1]],x)))==TRUE)==FALSE)
    stop("Individuals have different names across levels. Please check original matrices")
  
  # Are all variables to be analyzed numeric?
  for(i in 1:length(datalist)){
  if(all(apply(datalist[[i]][initialrow:nrow(datalist[[i]]),initialcolumn:ncol(datalist[[i]])],2,function(x) is.numeric(x))==F)) stop(paste("Check your input. Non numeric values in",names(datalist)[i], "dataset",sep=" "))
  }
  
  # Set imputation and abundance balancing defaults if not given by user
  if(is.null(imputation)) imputation <-"RF"
  if(is.null(abdbal)) abdbal <-"AvgIntensity"
  
  # Check imputation method(s)
  if(!is.null(imputation)) if((length(imputation)==1)&(imputation=="RF"|imputation=="KNN"|imputation=="none")) imputation <- rep(imputation, length(datalist))
  if(!is.null(imputation)) if((length(imputation)!=1)&(length(imputation)!=length(datalist))) stop("A vector containing one-common- or n elements (where n is the number of datasets) indicating imputation method(s) is expected")
  if(!is.null(imputation)) if(FALSE %in% (imputation %in% c("RF","KNN","none"))==TRUE) stop("Please select a valid imputation method")
  if(!is.null(imputation)) if(length(imputation)==1) imputation <- rep(imputation, length(datalist))
  
  # Check abundance balancing methods
  if(!is.null(abdbal)) if((length(abdbal)==1)&(abdbal=="sample"|abdbal=="AvgIntensity"|abdbal=="TreatAvgIntensity"|abdbal=="none")) abdbal <- rep(abdbal, length(datalist))
  if(!is.null(abdbal)) if(length(abdbal)!=1&length(abdbal)!=length(datalist)) stop("A vector containing one-common- or n elements (where n is the number of datasets) indicating abundance balancing method(s) is expected")
  if(!is.null(abdbal)) if(FALSE %in% (abdbal %in% c("sample","AvgIntensity","TreatAvgIntensity","none"))==TRUE) stop("Please select a valid abundance balancing method")
  if(!is.null(abdbal)) if(length(abdbal)==1) abdbal <- rep(abdbal, length(datalist))
  
  # Check k and threshold
  if(class(k)!="numeric") stop("Number of neighbors, k,  should be a numeric constant")
  #if(imputation=="RF") cat("Random Forest imputation. k value will be ignored") #Por que da esto error?
  
  if(class(threshold)!="numeric") stop("Threshold for defining NA or 0, threshold,  should be a numeric constant betwen 0 and 1")
  
  #### Imputation ####
  datasetnames<-names(datalist)

  if(parallel==FALSE){
    if(!is.null(imputation)){
      cat("\nMISSING VALUE IMPUTATION\n")
      cat("Single processor core will be used. It may take a while...\n")
      cat(paste(imputation, " imputation method will be used for ", names(datalist)," dataset\n",sep=""))
      datalist<-mapply(function(x,y,z){
        NAorZero_Imput(x,initialrow=initialrow,initialcolumn=initialcolumn,treatment1col=treatment1col, treatment2col=treatment2col, treatment=treatment,threshold=threshold,imputation=y,k=k,cores=1,datasetname = z)
      },datalist,imputation,datasetnames,SIMPLIFY = FALSE)
    }
  }

  if(parallel==TRUE){
    require(doParallel)
    cores=detectCores()
    if(!is.null(imputation)){
      cat("\nMISSING VALUE IMPUTATION\n")
      cat("Multiple processor cores (", cores, ") will be used. It may take a while...\n",sep="")
      cat(paste(imputation, " method will be used for ", names(datalist)," dataset\n",sep=""))
      datalist<-mapply(function(x,y,z){
        NAorZero_Imput(x,initialrow=initialrow,initialcolumn=initialcolumn,treatment1col=treatment1col, treatment2col=treatment2col, treatment=treatment,threshold=threshold,imputation=y,k=k,cores=cores,datasetname = z)
      },datalist,imputation,datasetnames,SIMPLIFY = FALSE)
    }
  }
  #### Removing empty columns ####
    cat(paste("\nREMOVING EMPTY COLUMNS OF ALL DATASETS\n"))
    cat("Single processor core will be used. It may take a while...\n")
    datalist<-lapply(datalist, function(x) RemoveEmptyColumns(x,initialrow=initialrow,initialcolumn=initialcolumn))
    names(datalist)<-datasetnames
  #### Abundance balancing ####
    if(!is.null(abdbal)){
      cat("\nABUNDANCE BALANCING\n")
      cat("Single processor core will be used. It may take a while...\n")
      cat(paste(abdbal, " balancing will be used for ", names(datalist)," dataset\n",sep=""))
      datalist<-mapply(function(x,y,z){
        AbdBal(x,initialrow=initialrow,initialcolumn=initialcolumn,treatment1col=treatment1col, treatment2col=treatment2col, treatment=treatment,norm=y,datasetname=z)
      },datalist,abdbal,datasetnames,SIMPLIFY = FALSE)
    }
  #### Final report ####
  header <- c("\n\n\nSUMMARY: Dataset imputation and balancing\n-----------------------------------------\n")
  header <- paste(header, paste(names(datalist),collapse=", "), " datasets were considered\n",sep="")

  missing <- c("\nMissing value imputation.........\n")
  if((imputation)==rep("none",length(imputation))){
    missing <- paste(missing,"Missing Value imputation was not performed.\n",sep="")
  } else if ("none" %in% imputation) {
    not<-which(imputation=="none")
    missing  <- paste(missing, paste(names(datalist)[-not],collapse=", "), " were respectively imputed employing ", paste(imputation[-not], collapse=", "), " methods.",sep="")
  } else{
    missing  <- paste(missing, paste(names(datalist),collapse=", "), " were respectively imputed employing ", paste(imputation, collapse=", "), " methods.",sep="")
    }
    
    if("KNN" %in% imputation) missing <- paste("\n",missing, k, "neighbors were employed for KNN calculations")
    missing <- paste(missing,"\nMinimum percentage of significant values per variable/treatment to allow imputation was ",threshold, sep=" ")
    missing <- paste(missing, "\nRun in parallel",parallel,"\n", sep=" ")
  
  if(is.null(abdbal)|(abdbal)==rep("none",length(abdbal))){
    abundance <- c("\nAbundance balancing.........\n")
    abundance <- paste(abundance, "Abundance balancing was not performed", sep="")
  } else if ("none" %in% abdbal) {
    abundance <- c("\nAbundance balancing.........\n")
    nein<-which(abdbal=="none")
    abundance <- paste(abundance, paste(names(datalist)[-nein], collapse=", "), " were respectively balanced employing ", paste(abdbal[-nein], collapse=", "), " methods.",sep="")
  } else{
    abundance <- c("\nAbundance balancing.........\n")
    abundance <- paste(abundance, paste(names(datalist), collapse=", "), " were respectively balanced employing ", paste(abdbal, collapse=", "), " methods.",sep="")
  }

  cat(header,missing,abundance)

  class(datalist) <-"POL"
  cat("\n\n Job finished!")
  return(datalist)

  
  }

RemoveEmptyColumns <-function(matriz, initialrow, initialcolumn){
  matrizfiltrada <-matriz[initialrow:nrow(matriz),initialcolumn:ncol(matriz)]
  vectorseleccion <- colSums(matrizfiltrada,na.rm=TRUE)==0
  if(length(vectorseleccion)==0){matrizfiltrada<-matrizfiltrada}
  else {matrizfiltrada <- matrizfiltrada[,vectorseleccion==F]}
  if(initialrow==1){
    matrizfiltrada <- cbind(matriz[,1:initialcolumn-1],matrizfiltrada)
    return(matrizfiltrada)
  } else if(initialrow>=2) {
    matrizfiltrada <- cbind(matriz[,1:initialcolumn-1],matrizfiltrada)
    matrizfiltrada <- rbind(matriz[1:initialrow-1,],matrizfiltrada)
    return(matrizfiltrada)
  } else {
    return(cat("\nError. Please revise your data input."))
  }
}

transformandselect <- function(datalist, initialrow=1, initialcolumn=3, treatment1col = 1, treatment2col = 2, treatment = 3, transf=NULL, varsel=FALSE,varselthld=0.25,varcoef=FALSE,varcoefthld=75){
  #### Initial checks ####
  if(class(datalist)!="POL") stop("Omic list should have a POL class. Please execute preprocess_omic_list first")
  
  # Set and check transformation 
  if(!is.null(transf)) if(length(transf)!=1&length(transf)!=length(datalist)) stop("A vector containing one-common- or n elements (where n is the number of datasets) indicating transformation method(s) is expected")
  if(!is.null(transf)) if(FALSE %in% (transf %in% c('z','z-centered','Log10', 'Log10plus1', 'sqrt', 'cubic', 'arcsine', 'box-cox-lambda','none'))==TRUE) stop("Please select a valid transformation method")
  if(!is.null(transf)) if(length(transf)==1) transf <- rep(transf, length(datalist))
  
  # Check variable selection settings
  if(varsel==TRUE & is.null(varselthld)) stop("Please select adequate threshold for VarSelect")
  if(!is.null(varselthld)) if(length(varselthld)!=1&length(varselthld)!=length(datalist)) stop("A vector containing one -common- or n (where n is the number of datasets) thresholds  is expected")
  if(!is.null(varselthld)) if(length(varselthld)==1) varselthld <- rep(varselthld, length(datalist))

  if(varcoef==TRUE & is.null(varcoefthld)) stop("Please select adequate threshold for VarCoefSelect")
  if(!is.null(varcoefthld)) if(length(varcoefthld)!=1&length(varcoefthld)!=length(datalist)) stop("A vector containing one -common- or n (where n is the number of datasets) thresholds  is expected")
  if(!is.null(varcoefthld)) if(length(varcoefthld)==1) varcoefthld <- rep(varcoefthld, length(datalist))

  datasetnames <-names(datalist)

  #### Transformation ####
  if(!is.null(transf)){
    cat("\nTRANSFORMATION OF DATASETS\n")
    cat("Single processor core will be used. It may take a while...\n")
    cat(paste(transf, " transformation method will be used for ", names(datalist)," dataset\n",sep=""))
    datalist<-mapply(function(x,y,z){
      Trnsfrm(x,initialrow=initialrow,initialcolumn=initialcolumn,treatment1col = treatment1col, treatment2col = treatment2col, treatment = treatment,datasetname = z,transform = y)
    },datalist,transf,datasetnames,SIMPLIFY = FALSE)
  }

  #### Variable selection ####
  if(varsel==TRUE){
    cat("\nSELECTING VARIABLES BASED ON CONSISTENCY\n")
    cat("Single processor core will be used. It may take a while...\n")
    datalist<-mapply(function(x,y,z){
      VarSelect(x,initialrow=initialrow,initialcolumn=initialcolumn,treatment1col = treatment1col, treatment2col = treatment2col, treatment = treatment,datasetname = z,threshold=y)
    },datalist,varselthld,datasetnames,SIMPLIFY = FALSE)
  }

  if(varcoef==TRUE){
    cat("\nSELECTING VARIABLES BASED ON COEFFICIENT OF VARIATION\n")
    cat("Single processor core will be used. It may take a while...\n")
    datalist<-mapply(function(x,y,z){
      VarCoefSelect(x,initialrow=initialrow,initialcolumn=initialcolumn,treatment1col = treatment1col, treatment2col = treatment2col, treatment = treatment,datasetname = z,threshold=y)
    },datalist,varcoefthld,datasetnames,SIMPLIFY = FALSE)
  }

  cat("\ndone!")

  class(datalist)<-"POL"

  return(datalist)
}

featuresel <-function(datalist, initialrow=1, initialcolumn=3, treatment1col=1, treatment2col=2, treatment=1, method=NULL, threshold=NULL, parametric=TRUE,stat="p"){
  if(class(datalist)!="POL") stop("Omic list should have a POL class. Please execute preprocess_omic_list first")

  # meter los checks correspondientes
  #if(!is.null(method)) if(method=="RF"|imputation=="KNN") imputation <- rep(imputation, length(datalist))

  #  if(!is.null(transf)) if(length(transf)!=1&length(transf)!=length(datalist)) stop("A vector containing one-common- or n elements (where n is the number of datasets) indicating transformation method(s) is expected")
  #  if(!is.null(transf)) if(FALSE %in% (transf %in% c('z','z-centered','Log10', 'log10plus1', 'sqrt', 'cubic', 'arcsine', 'box-cox-lambda'))==TRUE) stop("Please select a valid transformation method")
  #  if(!is.null(transf)) if(length(transf)==1) transf <- rep(transf, length(datalist))

  datasetnames <- names(datalist)
  if(!is.null(method)){
    cat("\nFEATURE SELECTION\n")
    cat("Single processor core will be used. It may take a while...\n")
    cat(paste(method, " selection method will be used for ", names(datalist)," dataset\n",sep=""))
    datalist<-mapply(function(x,y,z,a,b,c){
      feature_selection_list_helper(x,initialrow=initialrow,initialcolumn=initialcolumn,treatment1col=treatment1col,treatment2col=treatment2col,treatment=treatment,method=y,threshold=z,parametric=a,stat=b,datasetname=c)
    },datalist,method,threshold,parametric,stat,datasetnames,SIMPLIFY = FALSE)
  }
  class(datalist)<-"POL"
  return(datalist)
}

feature_selection_list_helper <- function(tabladatos, initialrow, initialcolumn, treatment1col, treatment2col, treatment, method, threshold, parametric,stat,datasetname=NULL){
  cat(paste("\nProcessing ",datasetname," dataset",sep=""))
  if(method=="IQR") tabladatos<-IQRfilter(tabladatos, initialrow = initialrow, initialcolumn=initialcolumn, treatment1col=treatment1col, treatment2col=treatment2col, treatment=treatment, method, threshold=threshold)
  if(method=="Stat") tabladatos<-Statfilter(tabladatos, initialrow = initialrow, initialcolumn=initialcolumn, treatment1col=treatment1col, treatment2col=treatment2col, initialcolumn=initialcolumn, threshold=threshold,parametric=parametric, stat=stat)
  return(tabladatos)
}

NAorZero_Imput <- function(matriz, initialrow, initialcolumn, treatment1col, treatment2col, treatment=1, threshold=0.25, imputation="RF",k=3, cores=2,datasetname=NULL){
  cat(paste("\nProcessing ",datasetname," dataset",sep=""))
  require("plyr")
  datos<-matriz[initialrow:nrow(matriz),initialcolumn:ncol(matriz)]
  initial0s <- sum(datos==0,na.rm=TRUE)
  initialNAs <-sum(is.na(datos))
  variablesinitial0s <- length(apply(apply(datos,2,function(x) x==0),2,any,na.rm=T)[apply(apply(datos,2,function(x) x==0),2,any,na.rm=T)==TRUE])
  variablesinitialNAs <- length(apply(apply(datos,2,function(x) is.na(x)),2,any)[apply(apply(datos,2,function(x) is.na(x)),2,any)==TRUE])
  datos[datos==0] <- NA
  matriz2 <-matriz
  matriz2[initialrow:nrow(matriz),initialcolumn:ncol(matriz)] <-datos #Reemplazamos por los datos con NA
  rm(datos) #vamos quitando basura
  #### Cero - NA Split ####
  if(treatment==1) vectortratamientos <- as.factor(matriz[initialrow:nrow(matriz),treatment1col])
  if(treatment==2) vectortratamientos <- as.factor(matriz[initialrow:nrow(matriz),treatment2col])
  if(treatment==3) vectortratamientos <- as.factor(paste(matriz[initialrow:nrow(matriz),treatment1col],matriz[initialrow:nrow(matriz),treatment2col], sep="&"))
  
  #ahora aqui hacemos un split, pero como llevamos todas las columnas anteriores, tampoco es muy critico
  datos <- split(as.data.frame(matriz2),vectortratamientos,drop=T) #Separo en tratamientos
  processedlist <- lapply(datos, function(x) NAorZeroTreatment(x,threshold,initialcolumn))
  matriz2 <- do.call("rbind",  unname(processedlist))
  #ahora reordenamos la matriz2 en base al nombre de filas de matriz1 para que se parezca a la original
  matriz2 <- matriz2[match(rownames(matriz), rownames(matriz2)), ]
  
  #ahora reordenamos la tabla para que se parezca a la inicial
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
    return(matriz2)
  } else if(imputation=="KNN"){
    require(impute)
    imputed<-matriz2[initialrow:nrow(matriz2),initialcolumn:ncol(matriz2)]
    imputed<-as.data.frame(impute.knn(as.matrix(imputed),k)$data)
    matriz2[initialrow:nrow(matriz2),initialcolumn:ncol(matriz2)]<-imputed
    
  } else if(imputation=="RF") {
    require(missForest)
    imputed<-matriz2[initialrow:nrow(matriz2),initialcolumn:ncol(matriz2)]
    if(cores>=2){
      require(doParallel)
      registerDoParallel(cores=cores) #require(doParallel)
      imputed<-(missForestedit(imputed,parallelize="variable"))$ximp
      matriz2[initialrow:nrow(matriz2),initialcolumn:ncol(matriz2)]<-imputed
    } else if(cores==1){
      imputed<-(missForest(imputed,parallelize="no"))$ximp}
    matriz2[initialrow:nrow(matriz2),initialcolumn:ncol(matriz2)]<-imputed
  } else{
    return(cat("\nError. Please revise selected imputation method."))
  }
  cat("Missing values have been imputed according to",imputation,"method.\n")
  return(matriz2)
}

NAorZeroTreatment <- function(submatriz,threshold,initialcolumn){
  submatrix <- submatriz[,initialcolumn:ncol(submatriz)]
  replicates <- nrow(submatrix)
  NACount <- apply(submatrix,2, function(x) sum(is.na(x)))
  Isa0 <- NACount <= (replicates*threshold)
  for(i in 1:length(Isa0)){
    if(Isa0[i]==FALSE){submatrix[,i][is.na(submatrix[,i])]<-0}
  }
  submatriz[,initialcolumn:ncol(submatriz)] <- submatrix
  return(submatriz)
}

missForestedit<-function (xmis, maxiter = 10, ntree = 100, variablewise = FALSE, decreasing = FALSE, verbose = FALSE, mtry = floor(sqrt(ncol(xmis))), replace = TRUE, classwt = NULL, cutoff = NULL, strata = NULL, 
                          sampsize = NULL, nodesize = NULL, maxnodes = NULL, xtrue = NA, parallelize = c("no", "variables", "forests")) {
  n <- nrow(xmis)
  p <- ncol(xmis)
  if (!is.null(classwt)) 
    stopifnot(length(classwt) == p, typeof(classwt) == "list")
  if (!is.null(cutoff)) 
    stopifnot(length(cutoff) == p, typeof(cutoff) == "list")
  if (!is.null(strata)) 
    stopifnot(length(strata) == p, typeof(strata) == "list")
  if (!is.null(nodesize)) 
    stopifnot(length(nodesize) == 2)
  if (any(apply(is.na(xmis), 2, sum) == n)) {
    indCmis <- which(apply(is.na(xmis), 2, sum) == n)
    xmis <- xmis[, -indCmis]
    p <- ncol(xmis)
    cat("  removed variable(s)", indCmis, "due to the missingness of all entries\n")
  }
  parallelize <- match.arg(parallelize)
  if (parallelize %in% c("variables", "forests")) {
    if (getDoParWorkers() == 1) {
      stop("You must register a 'foreach' parallel backend to run 'missForest' in parallel. Set 'parallelize' to 'no' to compute serially.")
    }
    else if (verbose) {
      if (parallelize == "variables") {
        cat("  parallelizing over the variables of the input data matrix 'xmis'\n")
      }
      else {
        cat("  parallelizing computation of the random forest model objects\n")
      }
    }
    if (getDoParWorkers() > p) {
      stop("The number of parallel cores should not exceed the number of variables (p=", 
           p, ")")
    }
  }
  ximp <- xmis
  xAttrib <- lapply(xmis, attributes)
  varType <- character(p)
  for (t.co in 1:p) {
    if (is.null(xAttrib[[t.co]])) {
      varType[t.co] <- "numeric"
      ximp[is.na(xmis[, t.co]), t.co] <- mean(xmis[, t.co], 
                                              na.rm = TRUE)
    }
    else {
      varType[t.co] <- "factor"
      max.level <- max(table(ximp[, t.co]))
      class.assign <- sample(names(which(max.level == summary(ximp[, 
                                                                   t.co]))), 1)
      if (class.assign != "NA's") {
        ximp[is.na(xmis[, t.co]), t.co] <- class.assign
      }
      else {
        while (class.assign == "NA's") {
          class.assign <- sample(names(which(max.level == 
                                               summary(ximp[, t.co]))), 1)
        }
        ximp[is.na(xmis[, t.co]), t.co] <- class.assign
      }
    }
  }
  NAloc <- is.na(xmis)
  noNAvar <- apply(NAloc, 2, sum)
  sort.j <- order(noNAvar)
  if (decreasing) 
    sort.j <- rev(sort.j)
  sort.noNAvar <- noNAvar[sort.j]
  nzsort.j <- sort.j[sort.noNAvar > 0]
  if (parallelize == "variables") {
    "%cols%" <- get("%dopar%")
    idxList <- as.list(isplitVector(nzsort.j, chunkSize = getDoParWorkers()))
  }
  Ximp <- vector("list", maxiter)
  iter <- 0
  k <- length(unique(varType))
  convNew <- rep(0, k)
  convOld <- rep(Inf, k)
  OOBerror <- numeric(p)
  names(OOBerror) <- varType
  if (k == 1) {
    if (unique(varType) == "numeric") {
      names(convNew) <- c("numeric")
    }
    else {
      names(convNew) <- c("factor")
    }
    convergence <- c()
    OOBerr <- numeric(1)
  }
  else {
    names(convNew) <- c("numeric", "factor")
    convergence <- matrix(NA, ncol = 2)
    OOBerr <- numeric(2)
  }
  stopCriterion <- function(varType, convNew, convOld, iter, 
                            maxiter) {
    k <- length(unique(varType))
    if (k == 1) {
      (convNew < convOld) & (iter < maxiter)
    }
    else {
      ((convNew[1] < convOld[1]) | (convNew[2] < convOld[2])) & 
        (iter < maxiter)
    }
  }
  while (stopCriterion(varType, convNew, convOld, iter, maxiter)) {
    if (iter != 0) {
      convOld <- convNew
      OOBerrOld <- OOBerr
    }
    cat("  missForest iteration", iter + 1, "in progress...")
    t.start <- proc.time()
    ximp.old <- ximp
    if (parallelize == "variables") {
      for (idx in idxList) {
        results <- foreach(varInd = idx, .packages = "randomForest") %cols% 
        {
          obsi <- !NAloc[, varInd]
          misi <- NAloc[, varInd]
          obsY <- ximp[obsi, varInd]
          obsX <- ximp[obsi, seq(1, p)[-varInd]]
          misX <- ximp[misi, seq(1, p)[-varInd]]
          typeY <- varType[varInd]
          if (typeY == "numeric") {
            set.seed(1511)
            RF <- randomForest(x = obsX, y = obsY, 
                               ntree = ntree, mtry = mtry, replace = replace, 
                               sampsize = if (!is.null(sampsize)) 
                                 sampsize[[varInd]]
                               else if (replace) 
                                 nrow(obsX)
                               else ceiling(0.632 * nrow(obsX)), nodesize = if (!is.null(nodesize)) 
                                 nodesize[1]
                               else 1, maxnodes = if (!is.null(maxnodes)) 
                                 maxnodes
                               else NULL)
            oerr <- RF$mse[ntree]
            misY <- predict(RF, misX)
          }
          else {
            obsY <- factor(obsY)
            summarY <- summary(obsY)
            if (length(summarY) == 1) {
              oerr <- 0
              misY <- factor(rep(names(summarY), length(misi)))
            }
            else {
              set.seed(1511)
              RF <- randomForest(x = obsX, y = obsY, 
                                 ntree = ntree, mtry = mtry, replace = replace, 
                                 classwt = if (!is.null(classwt)) 
                                   classwt[[varInd]]
                                 else rep(1, nlevels(obsY)), cutoff = if (!is.null(cutoff)) 
                                   cutoff[[varInd]]
                                 else rep(1/nlevels(obsY), nlevels(obsY)), 
                                 strata = if (!is.null(strata)) 
                                   strata[[varInd]]
                                 else obsY, sampsize = if (!is.null(sampsize)) 
                                   sampsize[[varInd]]
                                 else if (replace) 
                                   nrow(obsX)
                                 else ceiling(0.632 * nrow(obsX)), nodesize = if (!is.null(nodesize)) 
                                   nodesize[2]
                                 else 5, maxnodes = if (!is.null(maxnodes)) 
                                   maxnodes
                                 else NULL)
              oerr <- RF$err.rate[[ntree, 1]]
              misY <- predict(RF, misX)
            }
          }
          list(varInd = varInd, misY = misY, oerr = oerr)
        }
        for (res in results) {
          misi <- NAloc[, res$varInd]
          ximp[misi, res$varInd] <- res$misY
          OOBerror[res$varInd] <- res$oerr
        }
      }
    }
    else {
      for (s in 1:p) {
        varInd <- sort.j[s]
        if (noNAvar[[varInd]] != 0) {
          obsi <- !NAloc[, varInd]
          misi <- NAloc[, varInd]
          obsY <- ximp[obsi, varInd]
          obsX <- ximp[obsi, seq(1, p)[-varInd]]
          misX <- ximp[misi, seq(1, p)[-varInd]]
          typeY <- varType[varInd]
          if (typeY == "numeric") {
            if (parallelize == "forests") {
              xntree <- NULL
              set.seed(1511)
              RF <- foreach(xntree = idiv(ntree, chunks = getDoParWorkers()), 
                            .combine = "combine", .multicombine = TRUE, 
                            .packages = "randomForest") %dopar% {
                              randomForest(x = obsX, y = obsY, ntree = xntree, 
                                           mtry = mtry, replace = replace, sampsize = if (!is.null(sampsize)) 
                                             sampsize[[varInd]]
                                           else if (replace) 
                                             nrow(obsX)
                                           else ceiling(0.632 * nrow(obsX)), nodesize = if (!is.null(nodesize)) 
                                             nodesize[1]
                                           else 1, maxnodes = if (!is.null(maxnodes)) 
                                             maxnodes
                                           else NULL)
                            }
              OOBerror[varInd] <- mean((predict(RF) - 
                                          RF$y)^2, na.rm = TRUE)
            }
            else {
              set.seed(1511)        
              RF <- randomForest(x = obsX, y = obsY, 
                                 ntree = ntree, mtry = mtry, replace = replace, 
                                 sampsize = if (!is.null(sampsize)) 
                                   sampsize[[varInd]]
                                 else if (replace) 
                                   nrow(obsX)
                                 else ceiling(0.632 * nrow(obsX)), nodesize = if (!is.null(nodesize)) 
                                   nodesize[1]
                                 else 1, maxnodes = if (!is.null(maxnodes)) 
                                   maxnodes
                                 else NULL)
              OOBerror[varInd] <- RF$mse[ntree]
            }
            misY <- predict(RF, misX)
          }
          else {
            obsY <- factor(obsY)
            summarY <- summary(obsY)
            if (length(summarY) == 1) {
              misY <- factor(rep(names(summarY), sum(misi)))
            }
            else {
              if (parallelize == "forests") {
                set.seed(1511)
                RF <- foreach(xntree = idiv(ntree, chunks = getDoParWorkers()), 
                              .combine = "combine", .multicombine = TRUE, 
                              .packages = "randomForest") %dopar% 
                              {
                                set.seed(1511)
                                randomForest(x = obsX, y = obsY, 
                                             ntree = xntree, mtry = mtry, replace = replace, 
                                             classwt = if (!is.null(classwt)) 
                                               classwt[[varInd]]
                                             else rep(1, nlevels(obsY)), cutoff = if (!is.null(cutoff)) 
                                               cutoff[[varInd]]
                                             else rep(1/nlevels(obsY), nlevels(obsY)), 
                                             strata = if (!is.null(strata)) 
                                               strata[[varInd]]
                                             else obsY, sampsize = if (!is.null(sampsize)) 
                                               sampsize[[varInd]]
                                             else if (replace) 
                                               nrow(obsX)
                                             else ceiling(0.632 * nrow(obsX)), 
                                             nodesize = if (!is.null(nodesize)) 
                                               nodesize[2]
                                             else 5, maxnodes = if (!is.null(maxnodes)) 
                                               maxnodes
                                             else NULL)
                              }
                ne <- as.integer(predict(RF)) != as.integer(RF$y)
                ne <- ne[!is.na(ne)]
                OOBerror[varInd] <- sum(ne)/length(ne)
              }
              else {
                set.seed(1511)
                RF <- randomForest(x = obsX, y = obsY, 
                                   ntree = ntree, mtry = mtry, replace = replace, 
                                   classwt = if (!is.null(classwt)) 
                                     classwt[[varInd]]
                                   else rep(1, nlevels(obsY)), cutoff = if (!is.null(cutoff)) 
                                     cutoff[[varInd]]
                                   else rep(1/nlevels(obsY), nlevels(obsY)), 
                                   strata = if (!is.null(strata)) 
                                     strata[[varInd]]
                                   else obsY, sampsize = if (!is.null(sampsize)) 
                                     sampsize[[varInd]]
                                   else if (replace) 
                                     nrow(obsX)
                                   else ceiling(0.632 * nrow(obsX)), nodesize = if (!is.null(nodesize)) 
                                     nodesize[2]
                                   else 5, maxnodes = if (!is.null(maxnodes)) 
                                     maxnodes
                                   else NULL)
                OOBerror[varInd] <- RF$err.rate[[ntree, 
                                                 1]]
              }
              misY <- predict(RF, misX)
            }
          }
          ximp[misi, varInd] <- misY
        }
      }
    }
    cat("done!\n")
    iter <- iter + 1
    Ximp[[iter]] <- ximp
    t.co2 <- 1
    for (t.type in names(convNew)) {
      t.ind <- which(varType == t.type)
      if (t.type == "numeric") {
        convNew[t.co2] <- sum((ximp[, t.ind] - ximp.old[, 
                                                        t.ind])^2)/sum(ximp[, t.ind]^2)
      }
      else {
        dist <- sum(as.character(as.matrix(ximp[, t.ind])) != 
                      as.character(as.matrix(ximp.old[, t.ind])))
        convNew[t.co2] <- dist/(n * sum(varType == "factor"))
      }
      t.co2 <- t.co2 + 1
    }
    if (!variablewise) {
      NRMSE <- sqrt(mean(OOBerror[varType == "numeric"])/var(as.vector(as.matrix(xmis[, 
                                                                                      varType == "numeric"])), na.rm = TRUE))
      PFC <- mean(OOBerror[varType == "factor"])
      if (k == 1) {
        if (unique(varType) == "numeric") {
          OOBerr <- NRMSE
          names(OOBerr) <- "NRMSE"
        }
        else {
          OOBerr <- PFC
          names(OOBerr) <- "PFC"
        }
      }
      else {
        OOBerr <- c(NRMSE, PFC)
        names(OOBerr) <- c("NRMSE", "PFC")
      }
    }
    else {
      OOBerr <- OOBerror
      names(OOBerr)[varType == "numeric"] <- "MSE"
      names(OOBerr)[varType == "factor"] <- "PFC"
    }
    if (any(!is.na(xtrue))) {
      err <- suppressWarnings(mixError(ximp, xmis, xtrue))
    }
    if (verbose) {
      delta.start <- proc.time() - t.start
      if (any(!is.na(xtrue))) {
        cat("    error(s):", err, "\n")
      }
      cat("    estimated error(s):", OOBerr, "\n")
      cat("    difference(s):", convNew, "\n")
      cat("    time:", delta.start[3], "seconds\n\n")
    }
  }
  if (iter == maxiter) {
    if (any(is.na(xtrue))) {
      out <- list(ximp = Ximp[[iter]], OOBerror = OOBerr)
    }
    else {
      out <- list(ximp = Ximp[[iter]], OOBerror = OOBerr, 
                  error = err)
    }
  }
  else {
    if (any(is.na(xtrue))) {
      out <- list(ximp = Ximp[[iter - 1]], OOBerror = OOBerrOld)
    }
    else {
      out <- list(ximp = Ximp[[iter - 1]], OOBerror = OOBerrOld, 
                  error = suppressWarnings(mixError(Ximp[[iter - 
                                                            1]], xmis, xtrue)))
    }
  }
  class(out) <- "missForest"
  return(out)
}


AbdBal <- function(matriz, norm="AvgIntensity", initialrow=1, initialcolumn=3, treatment1col=1, 
                   treatment2col=2, treatment=1,datasetname=NULL){
  cat(paste("\nProcessing ",datasetname," dataset",sep=""))
  if(norm=="none"){ #Esto es no hacer nada.
    return(matriz)
  }
  #### Sample ####
  else if(norm=="sample"){ #Esto es lo mas facil.
    datos<-matriz[initialrow:nrow(matriz),initialcolumn:ncol(matriz)] #Saco los datos del resto
#    abundanciamedia <- mean(rowSums(datos,na.rm=T),na.rm=T)
    matriznormalizada <- t(apply(datos,1,function(x) x/sum(x, na.rm=T)))
    matriz[initialrow:nrow(matriz),initialcolumn:ncol(matriz)] <-matriznormalizada
    return(matriz)
  }
  #### AvgIntensity ####
  else if(norm=="AvgIntensity"){ #Esto es lo segundo mas facil. Centramos por cada muestra y multiplicamos por la media total
    datos<-matriz[initialrow:nrow(matriz),initialcolumn:ncol(matriz)] #Saco los datos del resto
    abundanciamedia <- mean(rowSums(datos,na.rm=T),na.rm=T)
    matriznormalizada <- t(apply(datos,1,function(x) x*abundanciamedia/sum(x, na.rm=T)))
    matriz[initialrow:nrow(matriz),initialcolumn:ncol(matriz)] <-matriznormalizada
    return(matriz)}
  #### TreatAvgIntensity ####
  else if(norm=="TreatAvgIntensity"){ 
    #Esto es lo tercero mas facil. Centramos por cada muestra y multiplicamos por la media total de cada tratamiento. Permitiria detectar tratamientos que causasen una disminucion de todas las variables
    #Esquema copiado de NA Split
    if(treatment==1) vectortratamientos <- as.factor(matriz[initialrow:nrow(matriz),treatment1col])
    if(treatment==2) vectortratamientos <- as.factor(matriz[initialrow:nrow(matriz),treatment2col])
    if(treatment==3) vectortratamientos <- as.factor(paste(matriz[initialrow:nrow(matriz),treatment1col],matriz[initialrow:nrow(matriz),treatment2col], sep="&"))
    datos <- split(as.data.frame(matriz),vectortratamientos,drop=T) #Separo en tratamientos
    processedlist <- lapply(datos, function(x) AbdBalTreatment(x,initialcolumn))
    matriz2 <- do.call("rbind",  unname(processedlist))
    #ahora reordenamos la matriz2 en base al nombre de filas de matriz1 para que se parezca a la original
    matriz2 <- matriz2[match(rownames(matriz), rownames(matriz2)), ]
    return(matriz2)
  }
}

AbdBalTreatment <- function(submatriz,initialcolumn){
  submatrix <- submatriz[,initialcolumn:ncol(submatriz)]
  abundanciamedia <- mean(rowSums(submatrix,na.rm=T),na.rm=T)
  submatrix <- t(apply(submatrix,1,function(x) x*abundanciamedia/sum(x, na.rm=T)))
  submatriz[,initialcolumn:ncol(submatriz)] <- submatrix
  return(submatriz)
}

Trnsfrm <- function(matriz, transformation="z", initialrow=1, initialcolumn=3,treatment1col=1, 
                    treatment2col=2, treatment=1,datasetname=NULL){
  cat(paste("\nProcessing ",datasetname," dataset",sep=""))
  datos<-matriz[initialrow:nrow(matriz),initialcolumn:ncol(matriz)]
  switch(transformation,
         "z"=datostransformados <- scale(datos,center=F),
         "z-centered"=datostransformados <- scale(datos),
         "Log10"={ceros<-datos==0
         datostransformados <- log(datos,10)
         datostransformados[ceros==T]<-0},
         "Log10plus1"=datostransformados<-log(datos+1,10),
         "box-cox-lambda"={
           if(treatment==1) vectortratamientos <- as.factor(matriz[initialrow:nrow(matriz),treatment1col])
           if(treatment==2) vectortratamientos <- as.factor(matriz[initialrow:nrow(matriz),treatment2col])
           if(treatment==3) vectortratamientos <- as.factor(paste(matriz[initialrow:nrow(matriz),treatment1col],matriz[initialrow:nrow(matriz),treatment2col], sep="&"))
           datostransformados<-apply(as.data.frame(datos),2,function(x) boxcoxvariable(x,vectortratamientos))
           },
         "sqrt"=datostransformados <- datos^(1/2),
         "cubic"=datostransformados <- datos^(1/3),
         "arcsine"=datostransformados <- asin(datos/100)*2/pi,
         "none"=datostransformados<-datos,
         return(cat("\nTransformation not recognized. Please use 'z','z-centered','Log10', 'log10plus1', 'sqrt', 'cubic', 'arcsine', 'box-cox-lambda' or 'none'.")))
  matriz[initialrow:nrow(matriz),initialcolumn:ncol(matriz)] <-datostransformados
  return(matriz)
}

boxcoxvariable <-function(vector,tratamiento){
  require("MASS")
  ceros <- vector==0 #sacamos las posiciones con ceros
  vector[vector==0] <- 0.000000000001
  lambda <- seq(-5,5,0.2)
  box <- boxcox(vector~as.factor(tratamiento),lambda,plotit=F) #se sacan los distintos modelos lambda en un rango amplio
  cox <- data.frame(box$x, box$y) #hacemos un data frame con los resultados
  cox2 <- cox[with(cox, order(-cox$box.y)),] #ordenamos de mayor a menor
  lambda <- cox2[1,"box.x"] #extraemos el lambda con la maximum likelihood logaritmica.
  if(lambda ==0){
    t_vector <-log(vector)
  } else{
    t_vector <- (vector^lambda)#-1/lambda   #(vector^lambda)-1/lambda Forma original de BoxCox
  }
  t_vector[ceros==T] <-0
  return(t_vector)
}

VarSelect <- function(matriz, initialrow, initialcolumn, treatment1col, treatment2col, treatment=1, threshold=0.25,datasetname=NULL){
  cat(paste("\nProcessing ",datasetname," dataset",sep=""))
  if(treatment==1) vectortratamientos <- as.factor(matriz[initialrow:nrow(matriz),treatment1col])
  if(treatment==2) vectortratamientos <- as.factor(matriz[initialrow:nrow(matriz),treatment2col])
  if(treatment==3) vectortratamientos <- as.factor(paste(matriz[initialrow:nrow(matriz),treatment1col],matriz[initialrow:nrow(matriz),treatment2col], sep="&"))
#### Consistentes ####
  # aqui saco vector con variables consistentes, me da igual de que tratamiento vengan
  # ojo que el vector consist empezara en la primera columna significativa y no en la primera de la tabla
  datos <- split(as.data.frame(matriz[initialrow:nrow(matriz),initialcolumn:ncol(matriz)]),vectortratamientos,drop=T) #Separo en tratamientos
  maxsignificantvalues <-lapply(datos,function(x) colSums(x !=0)/nrow(x)) #Saco el numero de valores significativos por variable y tratamiento y divido por el 
                                                                            #numero de filas (replicas). Si 1, esta presente en todos
  consist <- do.call(pmax,maxsignificantvalues) # me quedo con el valor maximo
  consist <- consist ==1 # True aquellos casos que son consistativos /esto pasa cuando reemplazas y no lees/
#### Cuantitativo ####
  # aqui saco vector de las variables q tengan n valores significativos (n>=threshold)
  # ojo que el vector cuantit empezara en la primera columna significativa y no en la primera de la tabla
  minreplicates <- round(nrow(matriz[initialrow:nrow(matriz),])*threshold,0) # Obtemenos el numero de casos y multiplicamos por threshold para obtener el numero minimo de replicas
  cuantit <- colSums(matriz[initialrow:nrow(matriz),initialcolumn:ncol(matriz)] !=0) #Sacamos valores no cero por columna (todos los tratamientos)
  cuantit <- cuantit >=minreplicates # True aquellos casos que son cuantitativos
#### Seleccion ####
  seleccion <- consist|cuantit # Juntamos todos los true
  cat("\nVariable Selection based on consistency. Variable was present at least on",minreplicates, " cases, or all replicates of a treatment" )
  cat("\nInitial Variables: ", ncol(matriz[,initialcolumn:ncol(matriz)]), " Selected Variables: ", length(seleccion[seleccion==T]), " Removed Variables: ", length(seleccion[seleccion==F]))
  matrizfiltrada <- cbind(as.data.frame(matriz[,1:initialcolumn-1]),as.data.frame(matriz[,initialcolumn:ncol(matriz)][,seleccion==T])) #Seleccionamos aquellas variables que queremos, conservando las columnas iniciales
  
  return(matrizfiltrada)
  
}

VarCoefSelect <- function(matriz, initialrow, initialcolumn, treatment1col, treatment2col, treatment=1, threshold=0.75,datasetname=NULL){
  options(warn=-1) #Si algun treatment es 0 en todos los casos, habra warning de inf, con esto no se muestra
  cat(paste("\nProcessing ",datasetname," dataset",sep=""))
  if(treatment==1) vectortratamientos <- as.factor(matriz[initialrow:nrow(matriz),treatment1col])
  if(treatment==2) vectortratamientos <- as.factor(matriz[initialrow:nrow(matriz),treatment2col])
  if(treatment==3) vectortratamientos <- as.factor(paste(matriz[initialrow:nrow(matriz),treatment1col],matriz[initialrow:nrow(matriz),treatment2col], sep="&"))
  
  datos <- split(as.data.frame(matriz),vectortratamientos,drop=T) #Separo en tratamientos
  coefvariacion <-lapply(datos,function(x) apply(x[initialrow:nrow(x),initialcolumn:ncol(x)],2,function(y) sd(y)/mean(y))) #Ojo, si un tratamiento es 0 tendra NA.
  coefvariacion <- do.call("rbind",coefvariacion) #transformamos la lista de vectores a tabla
  coefvariacion <- apply(coefvariacion,2,min,na.rm=T) #nos quedamos con el valor mas bajo. 
  #Lau:esto habra que especificarlo: pasa el filtro si los valores tienen una dispersion razonable (osea pequeña) en al menos en alguna treatment category
  #Lau:esto a mi me parece muy laxo pero vale
  seleccion <- coefvariacion <= threshold
  cat("\nVariable Selection based on coefficient of variation.")
  cat("\nOut of ",length(seleccion),"variables,",length(seleccion[seleccion==T])," show a coefficient of variation below ",threshold,"% in at least one treatment category.")
  
  matrizfiltrada <- cbind(matriz[,1:initialcolumn-1], matriz[,initialcolumn:ncol(matriz)][,seleccion]) #Seleccionamos aquellas variables que queremos, conservando las columnas iniciales
  
  return(matrizfiltrada)
}

IQRfilter <- function(matriz, initialrow, initialcolumn, treatment1col, treatment2col, treatment=1, threshold=0.25,method){
  # Aqui no se cuida el orden tras split porque es irrelevante
  if(treatment!=3&treatment!=0){
    if(treatment==1) treatmentcolumn <-treatment1col
    if(treatment==2) treatmentcolumn <-treatment2col
    if(unique(matriz[,treatmentcolumn])==1) stop("More than one level is required")
    datos <- split(as.data.frame(matriz),matriz[,treatmentcolumn],drop=T)
    mediaportratamiento <-lapply(datos,function(x) colMeans(data.matrix(x[initialrow:nrow(x),initialcolumn:ncol(x)]),na.rm=TRUE)) #Sacamos los valores medios de cada tratamiento
    mediaportratamiento <- do.call(rbind,mediaportratamiento) # hacemos una tabla con tantas filas como tratamientos con la media de cada variable

  }
  if(treatment==3){
    splitvector <- as.vector(paste(matriz[initialrow:nrow(matriz),treatment1col], matriz[initialrow:nrow(matriz),treatment2col],sep="&"))
    datos <- split(as.data.frame(matriz),splitvector,drop=T)
    mediaportratamiento <-lapply(datos,function(x) colMeans(data.matrix(x[initialrow:nrow(x),initialcolumn:ncol(x)]),na.rm=TRUE)) #Sacamos los valores medios de cada tratamiento
    mediaportratamiento <- do.call(rbind,mediaportratamiento) # hacemos una tabla con tantas filas como tratamientos
  }
  if(treatment==0){ # se toman valores absolutos
    mediaportratamiento <-data.matrix(matriz[initialrow:nrow(matriz),initialcolumn:ncol(matriz)])
  }
  
  iqr <- apply(mediaportratamiento,2,function(x) IQR(x,type = 8)/mean(x)) # sacamos el IQR y lo dividimos por la media para evitar sesgos de magnitud
  seleccion<- iqr > mean(iqr)*(1+threshold) #seleccionamos como T aquellos casos que cumplen el criterio
  cat("\nVariable Selection based on IQR.")
  cat("\nOut of ",length(seleccion),"variables,",length(seleccion[seleccion==T])," showed a IQR",threshold*100,"% greater than mean IQR.")
 
  matrizfiltrada <- cbind(matriz[,1:initialcolumn-1], matriz[,initialcolumn:ncol(matriz)][,seleccion]) #Seleccionamos aquellas variables que queremos, conservando las columnas iniciales
  
  return(matrizfiltrada)
}

Statfilter <- function(matriz, initialrow=1, initialcolumn=3, treatment1col=1, treatment2col=2, treatment=1, threshold=0.25, parametric=TRUE, stat="p"){
  datos<-matriz[initialrow:nrow(matriz),initialcolumn:ncol(matriz)] #Saco los datos del resto
  if(treatment==1) tratamientos <- as.factor(matriz[initialrow:nrow(matriz),treatment1col])
  if(treatment==2) tratamientos <- as.factor(matriz[initialrow:nrow(matriz),treatment2col])
  if(treatment==3) tratamientos <- as.factor(paste(matriz[initialrow:nrow(matriz),treatment1col],matriz[initialrow:nrow(matriz),treatment2col], sep="&"))

  if(parametric==TRUE){
    pvalores<-sapply(datos, function(x) anova(aov(x~tratamientos))$'Pr(>F)'[1])
    message<-"ANOVA"
  } else {
    pvalores<-sapply(datos, function(x) kruskal.test(x~tratamientos)$p.value)
    message<-"Kruskal"
  }
  if(stat=="p"){
    seleccion <- pvalores<=threshold
  } else if(stat=="q") {
    qvalores <-p.adjust(pvalores, method = "BH")
    seleccion <- qvalores<=threshold
  } else {
    return(cat("\nError. stat = c('p','q'), Please select mode p or q value"))
  }
  matrizfiltrada<-datos[,seleccion]
  cat("\nOut of",ncol(datos),"initial variables,",ncol(matrizfiltrada), "variables have been selected according to",message ,"test, applying a",stat,"value less than or equal to defined",threshold,"threshold.")
  if(initialrow==1){
    matrizfiltrada <- cbind(matriz[,1:initialcolumn-1],matrizfiltrada)
    return(matrizfiltrada)
  } else if(initialrow>=2) {
    matrizfiltrada <- cbind(matriz[,1:initialcolumn-1],matrizfiltrada)
    matrizfiltrada <- rbind(matriz[1:initialrow-1,],matrizfiltrada)
    return(matrizfiltrada)
  } else {
    return(cat("\nError. Please revise your data input."))
  }
}
