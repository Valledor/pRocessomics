#
# Se hace con las feautred proteins/metabolites
#
# Los objetos ya estan creados: Needles_Filtr_Z Needles_Filtr_featured_Z

#
# Empezamos solo con las proteinas, discriminante usando origenes como factor
#
da_analysis <- function(listadatos, annotation=NULL, initialrow=1, initialcolumn=3, treatment1col=1, treatment2col=2, 
                        treatment=1, omiclevel="all", keepX=NULL, autotune=FALSE, performance=FALSE, folds=3){
  options(warn=-1)
  require(mixOmics)
  require(plyr)
  require(doParallel)
  
  if(class(listadatos)!="POL") stop("A POL object is required. Please run pre-processing steps")
  datasetnames<-names(listadatos)
  
  #Generamos un vector para hacer el split por tratamientos
  if(treatment%in%c(1,2,3)==FALSE) stop("Select a valid treatment")
  if(treatment==1) vectortratamientos <- as.character.factor(listadatos[[1]][initialrow:nrow(listadatos[[1]]),treatment1col])
  if(treatment==2) vectortratamientos <- as.character.factor(listadatos[[1]][initialrow:nrow(listadatos[[1]]),treatment2col])
  if(treatment==3) vectortratamientos <- as.character.factor(paste(listadatos[[1]][initialrow:nrow(listadatos[[1]]),treatment1col],listadatos[[1]][initialrow:nrow(listadatos[[1]]),treatment2col], sep="&"))
  
  #aprovechamos para hacer la tabla de tratamientos. si no se hace esto los plots pierden versatilidad y o petan
  if(is.null(treatment2col)) treatment2col<-treatment1col
  if(treatment1col==treatment2col) tratamientos <-cbind(as.character.factor(listadatos[[1]][initialrow:nrow(listadatos[[1]]),treatment1col]),as.character.factor(listadatos[[1]][initialrow:nrow(listadatos[[1]]),treatment1col]))
  if(treatment1col!=treatment2col) tratamientos <-cbind(as.character.factor(listadatos[[1]][initialrow:nrow(listadatos[[1]]),treatment1col]),as.character.factor(listadatos[[1]][initialrow:nrow(listadatos[[1]]),treatment2col]))
  tratamientos <- cbind(tratamientos, paste(as.character.factor(listadatos[[1]][initialrow:nrow(listadatos[[1]]),treatment1col]),as.character.factor(listadatos[[1]][initialrow:nrow(listadatos[[1]]),treatment2col]), sep="&"))
  
  #Extraemos los valores significativos de las tablas
  listadatos2 <- lapply(listadatos, function(x) x[initialrow:nrow(x),initialcolumn:ncol(x)])
  if(is.null(omiclevel)==TRUE) omiclevel="all"
  
  # Seleccionamos los niveles omicos datos
  if(omiclevel=="all") {
    variablenames <- lapply(listadatos2, function(x) colnames(x))
    variablenames <-unlist(variablenames)
    listaintermedia <-do.call(cbind,listadatos2)  #### Esto es inconsistente. Cambia nombres de columnas, a veces, nombretabla.nombrevar, otras veces solo nombrevar
    colnames(listaintermedia) <- variablenames #### para arreglar lode arriba entre tanto.
    listaintermedia<-RemoveEmptyColumns(listaintermedia,1,1)#Esto es para evitar cosas raras.
  }  
  if(omiclevel!="all") {
    sapply(omiclevel, function(x) if(x %in% datasetnames==FALSE) stop("Error: Please select a valid datasetname"))
    
    if(length(omiclevel)>1){  
    listaintermedia <-do.call(cbind,listadatos2[omiclevel])
    variablenames <- lapply(listadatos2[omiclevel], function(x) colnames(x))
    variablenames <-unlist(variablenames)
    colnames(listaintermedia) <- variablenames
    listaintermedia<-RemoveEmptyColumns(listaintermedia,1,1)#Esto es para evitar cosas raras.
    }
    if(length(omiclevel)==1){
      listaintermedia <-listadatos2[[omiclevel]]
      variablenames <-colnames(listaintermedia)
      listaintermedia<-RemoveEmptyColumns(listaintermedia,1,1)#Esto es para evitar cosas raras.
    }
  }
  
  if(!is.null(annotation)){ #Creamos una tabla de anotaciones, a nuestra manera, para plotear luego. Chapuza, pero funciona.
    annotation_table <- annotation_mapman(listaintermedia,annotation)
  }
  
  
  if(autotune==TRUE){
  ###
  ### Facil de paralelizar pero lo dejo pendiente.
  ###
  miperf <- c()
  vectordepaso <-series(listaintermedia)

  pb <- txtProgressBar(min = 0, max = length(vectordepaso), style = 3)
  cat("\n Testing models for different KeppXs. These values will be tested:\n ",vectordepaso,"\n Single core will be used (by now) so it may take a while...")
  for(i in 1:length(vectordepaso)){
    set.seed(5881)
    da_tune <- splsda(listaintermedia,as.factor(vectortratamientos),mode ="regression",ncomp=2,keepX=c(vectordepaso[i],vectordepaso[i]))
    set.seed(5881)
    invisible(perform <- perf(da_tune, validation="Mfold", folds=3, criterion="Mfold", progressbar=FALSE,nrepeat = 50))
    miperf <-cbind(miperf, perform$error.rate$BER[,1])
    setTxtProgressBar(pb, i)
  }
  xc1<-min(which(miperf[1,]==min(miperf[1,]))) #seleccionamos el menor numero de variables para tener el error mas bajo en Comp1
  xc2<-min(which(miperf[2,]==min(miperf[2,]))) #seleccionamos el menor numero de variables para tener el error mas bajo en Comp2
  keepX <-c(xc1,xc2)
  cat("Best errors in components achieved for:\n", "Errors: ",miperf[keepX],"When selecting ",vectordepaso[keepX],  "Elements for Comp1 and Comp2 respectively")
  readline(prompt="Press [enter] to continue")
  }
  
  # Hacemos el análisis 
  cat("\nCALCULATING DA\nIt may take a while...\n")
  set.seed(5881)
  da_analysis <- splsda(listaintermedia,as.factor(vectortratamientos),mode ="regression",ncomp=2,keepX =keepX)
  
  if(performance==TRUE){
  cat("\nAnalyzing model performance\nIt may take a while...\n")
  tentfolds <- min(plyr::count(as.vector(vectortratamientos))[,2]) #Sacamos el fold maximo para que no pete performance
  if(is.null(folds)) folds<-tentfolds
  if(folds>tentfolds) folds<-tentfolds
  set.seed(5881)
  perform <- perf(da_analysis, validation="Mfold", folds=folds, criterion="Mfold", progressbar=TRUE,nrepeat = 50)
  plot(perform, col = color.mixo(5:7), sd = T, legend.position = "horizontal") 
  #if(askYesNo("Plot ROC plot?")==TRUE) auroc(da_analysis, roc.comp = 2) 
  }
  
  if(is.null(annotation)){
    resultado <- list(da_analysis,tratamientos,names(listadatos))
    names(resultado) <- c("da","treatments","datasetnames")
    class(resultado) <- "daanalysis"
    return(resultado)
  } else {
    #Esto es una chapuza para que las redes puedan salir en pantalla con nombres de variables
    colnames(listaintermedia) <- make.names(annotation_table[,2],unique = TRUE)
    set.seed(5881)
    da_network <- splsda(listaintermedia,as.data.frame(as.factor(vectortratamientos)),mode ="regression",ncomp=2,keepX =keepX)
    resultado <- list(da_analysis,tratamientos,names(listadatos),annotation_table,da_network)
    names(resultado) <- c("da","treatments","datasetnames","annotations","danetwork")
    class(resultado) <- "daanalysis"
    return(resultado)
  }
  
}

spls_analysis <- function(listadatos, annotation=NULL, initialrow=1, initialcolumn=3, treatment1col=1, treatment2col=2, 
                          treatment=1, omiclevelX,omiclevelY, keepX=NULL, keepY=NULL,autotune=FALSE, performance=FALSE, folds=3){
  options(warn=-1)
  require(mixOmics)
  require(plyr)
  
  if(class(listadatos)!="POL") stop("A POL object is required. Please run pre-processing steps")
  datasetnames<-names(listadatos)
  if(is.null(omiclevelX)) stop("Select omic level for prediction matrix (omiclevelX)")
  if(is.null(omiclevelY)) stop("Select omic level for prediction matrix (omiclevelY)")
  
  #Treatment table, needed for versatile plots
  if(is.null(treatment2col)) treatment2col<-treatment1col
  if(treatment1col==treatment2col) tratamientos <-cbind(listadatos[[1]][initialrow:nrow(listadatos[[1]]),treatment1col],listadatos[[1]][initialrow:nrow(listadatos[[1]]),treatment1col])
  if(treatment1col!=treatment2col) tratamientos <-cbind(listadatos[[1]][initialrow:nrow(listadatos[[1]]),treatment1col],listadatos[[1]][initialrow:nrow(listadatos[[1]]),treatment2col])
  tratamientos <- cbind(tratamientos, paste(listadatos[[1]][initialrow:nrow(listadatos[[1]]),treatment1col],listadatos[[1]][initialrow:nrow(listadatos[[1]]),treatment2col], sep="&"))
  
  #Extraction of numeric values to be analyzed (range defined by initial-end rows and columns)
  listadatos2 <- lapply(listadatos, function(x) x[initialrow:nrow(x),initialcolumn:ncol(x)])
  
  #### Building of predictor and response matrices ####
  sapply(omiclevelX, function(x) if(x %in% datasetnames==FALSE) stop("Error: Please select a valid datasetname for prediction matrix"))
  sapply(omiclevelY, function(x) if(x %in% datasetnames==FALSE) stop("Error: Please select a valid datasetname for response matrix"))
  
  if(length(omiclevelX)>1){  
    predictormatrix <-do.call(cbind,listadatos2[omiclevelX])
    colnames(predictormatrix) <- unlist(lapply(listadatos2[omiclevelX],function(x) colnames(x)))
    predictormatrix<-RemoveEmptyColumns(predictormatrix,1,1)#Esto es para evitar cosas raras.
  }
  if(length(omiclevelX)==1){
    predictormatrix <-listadatos2[[omiclevelX]]
    predictormatrix<-RemoveEmptyColumns(predictormatrix,1,1)#Esto es para evitar cosas raras.
  }
  
  if(length(omiclevelY)>1){  
    responsematrix <-do.call(cbind,listadatos2[omiclevelY])
    colnames(responsematrix) <- unlist(lapply(listadatos2[omiclevelY],function(x) colnames(x)))
    responsematrix<-RemoveEmptyColumns(responsematrix,1,1)#Esto es para evitar cosas raras.
  }
  if(length(omiclevelY)==1){
    responsematrix <-listadatos2[[omiclevelY]]
    responsematrix<-RemoveEmptyColumns(responsematrix,1,1)#Esto es para evitar cosas raras.
  }
  
  #### Building annotation table ####
  if(!is.null(annotation)){ #Creamos una tabla de anotaciones, a nuestra manera, para plotear luego. Chapuza, pero funciona.
    annotation_table_X <- annotation_mapman(predictormatrix,annotation)
    annotation_table_Y <- annotation_mapman(responsematrix,annotation)
    annotation_table <- rbind(annotation_table_X,annotation_table_Y)
    variablenamesX <- make.names(annotation_table_X[,2],unique=TRUE)
    variablenamesY <- make.names(annotation_table_Y[,2],unique=TRUE)
  }
  
  #### Autotune SPLS ####
  if(autotune==TRUE){
    ###
    ### Facil de paralelizar pero lo dejo pendiente. En este caso habría que probar distintos X e Y para los dos componentes
    ###
    miperf <- c()
    ##aqui habria que mejorar las reglas para hacer la secuencia, pero de momento funciona
    vectordepasoX <-series(predictormatrix)
    vectordepasoY <-series(responsematrix)
    
    pb <- txtProgressBar(min = 0, max = length(vectordepasoX), style = 3)
    cat("Testing models for different keepX in prediction and response matrices:\n ","\n Single core will be used (by now) so it may take a while... go for coffee")
    
    for(i in 1:length(vectordepasoX)){
      for(j in 1:length(vectordepasoY)){
        set.seed(5881)
        spls_tune <- spls(predictormatrix,responsematrix,mode ="regression",ncomp=2,keepX=c(vectordepasoX[i],vectordepasoX[i]),keepY = c(vectordepasoY[j],vectordepasoY[j]))
        set.seed(5881)
        perf <- perf(spls_tune, validation="Mfold", folds=3, criterion="Mfold", progressbar=FALSE,nrepeat = 50)
        miperf <- c(miperf,sum(perf$Q2.total))
      }
      setTxtProgressBar(pb, i)
    }
    mimatriz <- matrix(miperf,ncol = length(vectordepasoY),byrow = T)
    row.names(mimatriz)<-vectordepasoX
    colnames(mimatriz)<-vectordepasoY
    
    selection<- as.vector(which(mimatriz == min(mimatriz), arr.ind = TRUE))
    selection<-c(vectordepasoX[selection[1]],vectordepasoY[selection[2]])
    keepX <- c(selection[1],selection[1])
    keepY <- c(selection[2],selection[2])
    
    cat("\nLowest error (total Q2):",min(miperf), "Was obtained when selecting: ", selection[1], " and ", selection[2],  "variables of prediction and response matrix respectively")
    readline(prompt="Press [enter] to continue")
    
  }
  
  #### sPLS analysis ####
  cat("\nCALCULATING sPLS\nIt may take a while...\n")
  set.seed(5881)
  spls_analysis <- spls(predictormatrix,responsematrix,mode ="regression",ncomp=2,keepX=keepX,keepY =keepY)
  
  #### Performance SPLS ####
  if(performance==TRUE){
    #Generamos un vector para determinar los folds
    if(treatment%in%c(1,2,3)==FALSE) stop("Select a valid treatment")
    if(treatment==1) vectortratamientos <- as.factor(listadatos[[1]][initialrow:nrow(listadatos[[1]]),treatment1col])
    if(treatment==2) vectortratamientos <- as.factor(listadatos[[1]][initialrow:nrow(listadatos[[1]]),treatment2col])
    if(treatment==3) vectortratamientos <- as.factor(paste(listadatos[[1]][initialrow:nrow(listadatos[[1]]),treatment1col],listadatos[[1]][initialrow:nrow(listadatos[[1]]),treatment2col], sep="&"))
    cat("\nAnalyzing model performance\nIt may take a while...\n")
    tentfolds <- min(plyr::count(as.vector(vectortratamientos))[,2]) #Sacamos el fold maximo para que no pete performance
    if(is.null(folds)) folds<-tentfolds
    if(folds>tentfolds) folds<-tentfolds
    set.seed(5881)
    perform <- perf(spls_analysis, validation="Mfold", folds=folds, criterion="Mfold", progressbar=TRUE,nrepeat = 50)
    cat("\nQ2 Total errors:\n")
    cat(perform$Q2.total)
  }
  
  #### Return results ####
  # without annotations
  if(is.null(annotation)){
    resultado <- list(spls_analysis,tratamientos,names(listadatos))
    names(resultado) <- c("splso","treatments","datasetnames")
    class(resultado) <- "splsanalysis"
    return(resultado)
  } else {
    # Annotations. This is horrible, but makes the trick to get on screen networks with names
    colnames(predictormatrix) <-variablenamesX
    colnames(responsematrix) <- variablenamesY
    set.seed(5881)
    spls_network <- spls(predictormatrix,responsematrix,mode ="regression",ncomp=2,keepX=keepX,keepY =keepY)
    resultado <- list(spls_analysis,tratamientos,names(listadatos),annotation_table,spls_network)
    names(resultado) <- c("splso","treatments","datasetnames","annotations","splsnetwork")
    class(resultado) <- "splsanalysis"
    return(resultado)
  }
  
}

da_plot <- function(da_list,treatment=1, plottype="Distance",variablespace="x",useannot=FALSE,fontsizes=c(14,10,16,12),cutoff=0.8){ 
  library("plotly") 
  ## Plottype: Aqui estan, o estaran, todos los graficos que pueden ser utiles. Hay unos cuantos, la gente lega no los usara, 
  ## pero es interesante mostrarlos.
  #             "Composite" <- Collage con los 2 plots mas representativos del analisis
  #             "Composite_Annotated" <- lo mismo pero con anotaciones sin leyenda. bueno para explorar datos
  #             "Distance" <- representa los individuos (1 punto por muestra), considerando todos los niveles. Unico que permitira 3D
  #             "Var" <- Muestrea las correlaciones de variables, permite escoger dataset y anotacion
  #             "Network" <- Network, a implementar
  #             "All" <- Devuelve una lista con todos los graficos anteriores.
  
  
  ## ********Esto hay que mejorarlo un poco controlando todos los niveles del objeto y el tipo de grafico*******
  ## Validamos que el objeto de verdad sea del tipo que tiene que ser:
  ##
  if(class(da_list)!= "daanalysis")
    stop("DA analysis object is expected, please run da_analysis first")
  if(plottype %in% c("Composite","Distance","Var","Network","All","Composite_Annotated")==FALSE)
    stop("Please select a valid plot type")
  if(variablespace%in%c("x","y","xy")==FALSE) stop("Select an adequate variable representation space. x, y, or xy ")
  if(useannot==T&&is.null(da_list$annotations)){
    cat(paste("\nAnnotation matrix not provided. Considering useannot=F\n", "Please check/repeat da_analysis() including an annotation matrix\n","Read help or follow the tutorial for more information.\n",sep=""))
    useannot=FALSE
  }
  #...
  
  ##
  ## Definimos los objetos necesarios para dibujar
  ##
  
  da_object <- da_list$da
  treatmentstable <- da_list$treatments
  annotations <- da_list$annotations
  datasetnames <- da_list$datasetnames
  treatmentnames <- colnames(treatmentstable)
  
  #Plot sencillo, Distancias
  if(plottype=="Distance"){
    if(variablespace=="x"){
      explainedvariance <-as.vector(da_object$explained_variance$X)
      distanceplotdata <- da_object$variates$X
      xaxislabel <- paste("Component 1-",round(explainedvariance[1]*100,2),"% of variance", sep=" ")
      yaxislabel <- paste("Component 2-",round(explainedvariance[2]*100,2),"% of variance", sep=" ")
    }
    if(variablespace=="y"){
      explainedvariance <-da_object$explained_variance$Y
      distanceplotdata <- da_object$variates$Y
      xaxislabel <- paste("Component 1-",round(explainedvariance[1]*100,2),"% of variance", sep=" ")
      yaxislabel <- paste("Component 2-",round(explainedvariance[2]*100,2),"% of variance", sep=" ")
    }
    if(variablespace=="xy"){
      distanceplotdata <- (da_object$variates$X + da_object$variates$Y)/2
      xaxislabel <- "XY-Variate 1"
      yaxislabel <- "XY-Variate 2"
    }
    
    if(treatment!=4){  #Asi seleccionamos todas las visualizaciones que solo tengan un factor a mostrar.
      if(treatment==1) tratamientocolores <-as.vector(treatmentstable[,1])
      if(treatment==2) tratamientocolores <-as.vector(treatmentstable[,2])
      if(treatment==3) tratamientocolores <-as.vector(treatmentstable[,3])
      
        da_distance_plot <- plot_ly(as.data.frame(distanceplotdata), 
                                  x = ~distanceplotdata[,1], 
                                  y = ~distanceplotdata[,2],
                                  hoverinfo ="text",
                                  text = paste("Sample:", rownames(distanceplotdata), sep=" "),
                                  showlegend = TRUE,
                                  type="scatter",
                                  mode = "markers", 
                                  color = factor(tratamientocolores,labels = c(unique(tratamientocolores))),
                                  colors = paleta_tratamientos[1:length(unique(tratamientocolores))],
                                  #symbol = factor(rep(c(1:nrow(mcoin_object$mcoa$lambda)),each=nrow(mcoin_object$mcoa$SynVar)), labels=datasetnames),
                                  symbols = simbolos[1:length(unique(tratamientocolores))], 
                                  marker = list(size = 11)) %>%
          layout(xaxis = list(title = xaxislabel,titlefont=list(size=fontsizes[1]),tickfont=list(size=fontsizes[2])),
                 yaxis = list(title = yaxislabel,titlefont=list(size=fontsizes[1]),tickfont=list(size=fontsizes[2])),
                 title = "DA Maximum Distance Plot",
                 titlefont = list(size=fontsizes[3]),
                 legend = list(font = list(size=fontsizes[4])))
        
        return(da_distance_plot)
      }
      
    if(treatment==4){  #Asi seleccionamos las visualizaciones con 2 factores a mostrar
      
        da_distance_plot <- plot_ly(as.data.frame(distanceplotdata), 
                                  x = ~distanceplotdata[,1], 
                                  y = ~distanceplotdata[,2],
                                  hoverinfo ="text",
                                  text = paste("Sample:", rownames(distanceplotdata), sep=" "),
                                  showlegend = TRUE,
                                  type="scatter",
                                  mode = "markers", 
                                  color = factor(treatmentstable[,1], labels = unique(treatmentstable[,1])),
                                  colors = paleta_tratamientos[1:length(unique(treatmentstable[,1]))],
                                  symbol = factor(treatmentstable[,2], labels=unique(treatmentstable[,2])),
                                  symbols = simbolos[1:length(unique(treatmentstable[,2]))], 
                                  marker = list(size = 11)) %>%
          layout(xaxis = list(title = xaxislabel,titlefont=list(size=fontsizes[1]),tickfont=list(size=fontsizes[2])),
                 yaxis = list(title = yaxislabel,titlefont=list(size=fontsizes[1]),tickfont=list(size=fontsizes[2])),
                 title = "DA Maximum Distance Plot",
                 titlefont = list(size=fontsizes[3]),
                 legend = list(font = list(size=fontsizes[4])))
        return(da_distance_plot)
      

        
      }
    }

  if(plottype=="Var"){
   vartable <-invisible(plotVar(da_object,plot=F)) #esto es una chapuza. habria que diseñar funcion para sacar esto. de momento uso mixomics
   vartable<-vartable[,1:3]
   if(variablespace=="x"){
     explainedvariance <-as.vector(da_object$explained_variance$X)
     xaxislabel <- paste("Component 1-",round(explainedvariance[1]*100,2),"% of variance", sep=" ")
     yaxislabel <- paste("Component 2-",round(explainedvariance[2]*100,2),"% of variance", sep=" ")
   }
   if(variablespace=="y"){
     explainedvariance <-da_object$explained_variance$Y
     xaxislabel <- paste("Component 1-",round(explainedvariance[1]*100,2),"% of variance", sep=" ")
     yaxislabel <- paste("Component 2-",round(explainedvariance[2]*100,2),"% of variance", sep=" ")
   }
   if(variablespace=="xy"){
     xaxislabel <- "XY-Variate 1"
     yaxislabel <- "XY-Variate 2"
   }
   
   if(useannot==FALSE){
      da_VarPlot <- plot_ly(as.data.frame(vartable), 
                             x = ~vartable[,1], 
                             y = ~vartable[,2],
                             text = paste("ID:", rownames(vartable)), 
                             type="scatter",
                             mode = "markers", 
                             marker = list(size = 6, opacity=1), 
                             colors =paleta_niveles[1:length(datasetnames)])
      
      da_VarPlot <- layout(da_VarPlot, xaxis = list(title = xaxislabel,titlefont=list(size=fontsizes[1]),tickfont=list(size=fontsizes[2])),
                                  yaxis = list(title = yaxislabel,titlefont=list(size=fontsizes[1]),tickfont=list(size=fontsizes[2])),
                                  title = "DA Variable Plot",
                                  titlefont = list(size=fontsizes[3]),
                                  legend = list(font = list(size=fontsizes[4])))
      
      
      return(da_VarPlot)
    }
    if(useannot==TRUE){
      #Generamos anotaciones
      
      annotations2 <- as.data.frame(row.names(vartable))
      colnames(annotations2)<-"IDENTIFIER"
      annotations <- left_join(annotations2,as.data.frame(annotations),by="IDENTIFIER") #con esto asignamos cada anotación a la variable correspondiente por rowname. Asi no hay fallo.
      annotations <- as.matrix(annotations) #Por alguna razon, sin esto, luego no se ordenan los factores ni los colores.
      rm(annotations2)
      
      #Sacamos los niveles del mapman, y las etiquetas
      mapmanlevels <- as.numeric(unique(annotations[order(as.numeric(annotations[,3])),3]))
      mapmannames <- annotations[order(as.numeric(annotations[,3])),]
      mapmannames <- as.vector(unique(mapmannames[,4]))
      
      da_VarPlot <- plot_ly(as.data.frame(vartable), 
                             x = ~vartable[,1], 
                             y = ~vartable[,2],
                             text = paste("ID:", rownames(vartable),"</br> Description: ",annotations[,2],"</br> Mercator Bin: ",annotations[,4]), 
                             type="scatter",
                             mode = "markers", 
                             marker = list(size = 6, opacity=1),
                             color = factor(annotations[,4], levels=mapmannames),
                             colors = paleta_mapman[mapmanlevels])
      
      da_VarPlot <- layout(da_VarPlot, xaxis = list(title = xaxislabel,titlefont=list(size=fontsizes[1]),tickfont=list(size=fontsizes[2])),
                           yaxis = list(title = yaxislabel,titlefont=list(size=fontsizes[1]),tickfont=list(size=fontsizes[2])),
                           title = "DA Variable Plot",
                           titlefont = list(size=fontsizes[3]),
                           legend = list(font = list(size=fontsizes[4])))
      
      return(da_VarPlot)
    }
  }
  
  if(plottype=="Composite"){
    p1 <- da_plot(da_list,treatment, plottype = "Distance", fontsizes = fontsizes)
    p2 <- da_plot(da_list,treatment, plottype = "Var", fontsizes = fontsizes)
    da_CompositeplotPlot <- subplot(subplot(p1,p2),nrows=2,shareX = F,shareY = F,margin = 0.04)
    da_CompositeplotPlot <- layout(da_CompositeplotPlot, title = "DA Analysis",titlefont = list(size=fontsizes[3]))
    return(da_CompositeplotPlot)
  }
  
  if(plottype=="Composite_Annotated"){
    p1 <- da_plot(da_list,treatment, plottype = "Distance", fontsizes = fontsizes)
    p2 <- da_plot(da_list,treatment, plottype = "Var", fontsizes = fontsizes,useannot = T)
    spls_CompositeplotPlot <- subplot(subplot(p1,p2),nrows=2,shareX = F,shareY = F,margin = 0.04)
    spls_CompositeplotPlot <- layout(spls_CompositeplotPlot, title = "SPLS Analysis",titlefont = list(size=fontsizes[3]),showlegend=FALSE)
    return(spls_CompositeplotPlot)
  }
  
  if(plottype=="Network"){
    if(useannot==FALSE) {
      set.seed(5881)
      network(da_object, comp=c(1,2),cutoff = cutoff)
      class(da_object) <- c(class(da_object),"networkplot")
      return(da_object)
      }
    if(useannot==TRUE){
      set.seed(5881)
      network(da_list[[5]], comp=c(1,2),cex.node.name = 0.5,cutoff = cutoff)
      class(da_list[[5]]) <- c(class(da_list[[5]]),"networkplotannot")
      return(da_list[[5]])
  }
  }
}

spls_plot <- function(spls_list,treatment=1, plottype="Distance",variablespace="xy",useannot=FALSE,fontsizes=c(14,10,16,12),cutoff=0.8){ 
  library("plotly") 
  ## Plottype: Aqui estan, o estaran, todos los graficos que pueden ser utiles. Hay unos cuantos, la gente lega no los usara, 
  ## pero es interesante mostrarlos.
  #             "Composite" <- Collage con los 2 plots mas representativos del analisis
  #             "Distance" <- representa los individuos (1 punto por muestra), considerando todos los niveles. Unico que permitira 3D
  #             "Var" <- Muestrea las correlaciones de variables, permite escoger dataset y anotacion
  #             "Network" <- Network, a implementar
  #             "All" <- Devuelve una lista con todos los graficos anteriores.
  
  
  ## ********Esto hay que mejorarlo un poco controlando todos los niveles del objeto y el tipo de grafico*******
  ## Validamos que el objeto de verdad sea del tipo que tiene que ser:
  ##
  if(class(spls_list)!= "splsanalysis")
    stop("SPLS analysis object is expected, please run spls_analysis first")
  if(plottype %in% c("Composite","Distance","DistanceXY", "Var","Network","All","Composite_Annotated")==FALSE)
    stop("Please select a valid plot type")
  if(variablespace%in%c("x","y","xy")==FALSE) stop("Select an adequate variable representation space. x, y, or xy ")
  if(useannot==T&&is.null(spls_list$annotations)){
    cat(paste("\nAnnotation matrix not provided. Considering useannot=F\n", "Please check/repeat da_analysis() including an annotation matrix\n","Read help or follow the tutorial for more information.\n",sep=""))
    useannot=FALSE
  }
  #...
  
  ##
  ## Definimos los objetos necesarios para dibujar
  ##
  
  spls_object <- spls_list$splso
  treatmentstable <- spls_list$treatments
  annotations <- spls_list$annotations
  datasetnames <- spls_list$datasetnames
  treatmentnames <- colnames(treatmentstable)
  
  #Plot sencillo, Distancias
  if(plottype=="Distance"){
    if(variablespace=="x"){
      explainedvariance <-as.vector(spls_object$explained_variance$X)
      distanceplotdata <- spls_object$variates$X
      xaxislabel <- paste("Component 1-",round(explainedvariance[1]*100,2),"% of variance", sep=" ")
      yaxislabel <- paste("Component 2-",round(explainedvariance[2]*100,2),"% of variance", sep=" ")
    }
    if(variablespace=="y"){
      explainedvariance <-spls_object$explained_variance$Y
      distanceplotdata <- spls_object$variates$Y
      xaxislabel <- paste("Component 1-",round(explainedvariance[1]*100,2),"% of variance", sep=" ")
      yaxislabel <- paste("Component 2-",round(explainedvariance[2]*100,2),"% of variance", sep=" ")
    }
    if(variablespace=="xy"){
      distanceplotdata <- (spls_object$variates$X + spls_object$variates$Y)/2
      xaxislabel <- "XY-Variate 1"
      yaxislabel <- "XY-Variate 2"
    }
    
    if(treatment!=4){  #Asi seleccionamos todas las visualizaciones que solo tengan un factor a mostrar.
      if(treatment==1) tratamientocolores <-as.vector(treatmentstable[,1])
      if(treatment==2) tratamientocolores <-as.vector(treatmentstable[,2])
      if(treatment==3) tratamientocolores <-as.vector(treatmentstable[,3])
      
      da_distance_plot <- plot_ly(as.data.frame(distanceplotdata), 
                                  x = ~distanceplotdata[,1], 
                                  y = ~distanceplotdata[,2],
                                  hoverinfo ="text",
                                  text = paste("Sample:", rownames(distanceplotdata), sep=" "),
                                  showlegend = TRUE,
                                  type="scatter",
                                  mode = "markers", 
                                  color = factor(tratamientocolores,labels = c(unique(tratamientocolores))),
                                  colors = paleta_tratamientos[1:length(unique(tratamientocolores))],
                                  #symbol = factor(rep(c(1:nrow(mcoin_object$mcoa$lambda)),each=nrow(mcoin_object$mcoa$SynVar)), labels=datasetnames),
                                  symbols = simbolos[1:length(unique(tratamientocolores))], 
                                  marker = list(size = 11)) %>%
        layout(xaxis = list(title = xaxislabel,titlefont=list(size=fontsizes[1]),tickfont=list(size=fontsizes[2])),
               yaxis = list(title = yaxislabel,titlefont=list(size=fontsizes[1]),tickfont=list(size=fontsizes[2])),
               title = "SPLS Maximum Distance Plot",
               titlefont = list(size=fontsizes[3]),
               legend = list(font = list(size=fontsizes[4])))
      
      return(da_distance_plot)
    }
    
    if(treatment==4){  #Asi seleccionamos las visualizaciones con 2 factores a mostrar
      
      da_distance_plot <- plot_ly(as.data.frame(distanceplotdata), 
                                  x = ~distanceplotdata[,1], 
                                  y = ~distanceplotdata[,2],
                                  hoverinfo ="text",
                                  text = paste("Sample:", rownames(distanceplotdata), sep=" "),
                                  showlegend = TRUE,
                                  type="scatter",
                                  mode = "markers", 
                                  color = factor(treatmentstable[,1], labels = unique(treatmentstable[,1])),
                                  colors = paleta_tratamientos[1:length(unique(treatmentstable[,1]))],
                                  symbol = factor(treatmentstable[,2], labels=unique(treatmentstable[,2])),
                                  symbols = simbolos[1:length(unique(treatmentstable[,2]))], 
                                  marker = list(size = 11)) %>%
        layout(xaxis = list(title = xaxislabel,titlefont=list(size=fontsizes[1]),tickfont=list(size=fontsizes[2])),
               yaxis = list(title = yaxislabel,titlefont=list(size=fontsizes[1]),tickfont=list(size=fontsizes[2])),
               title = "SPLS Maximum Distance Plot",
               titlefont = list(size=fontsizes[3]),
               legend = list(font = list(size=fontsizes[4])))
      return(da_distance_plot)
    }
  }
  
  if(plottype=="Var"){
    vartable <-invisible(plotVar(spls_object,plot=FALSE)) #esto es una chapuza. habria que diseñar funcion para sacar esto. de momento uso mixomics
    vartable<-vartable[,1:3]
    if(variablespace=="x"){
      explainedvariance <-as.vector(spls_object$explained_variance$X)
      xaxislabel <- paste("Component 1-",round(explainedvariance[1]*100,2),"% of variance", sep=" ")
      yaxislabel <- paste("Component 2-",round(explainedvariance[2]*100,2),"% of variance", sep=" ")
    }
    if(variablespace=="y"){
      explainedvariance <-spls_object$explained_variance$Y
      xaxislabel <- paste("Component 1-",round(explainedvariance[1]*100,2),"% of variance", sep=" ")
      yaxislabel <- paste("Component 2-",round(explainedvariance[2]*100,2),"% of variance", sep=" ")
    }
    if(variablespace=="xy"){
      xaxislabel <- "XY-Variate 1"
      yaxislabel <- "XY-Variate 2"
    }
    
    if(useannot==FALSE){
      spls_VarPlot <- plot_ly(as.data.frame(vartable), 
                              x = ~vartable[,1], 
                              y = ~vartable[,2],
                              text = paste("ID:", rownames(vartable)), 
                              type="scatter",
                              mode = "markers", 
                              marker = list(size = 6, opacity=1), 
                              colors =paleta_niveles[1:length(datasetnames)],
                              symbol = factor(vartable[,3], labels=unique(vartable[,3])),
                              symbols = simbolos[1:length(unique(vartable[,3]))])
      
      spls_VarPlot <- layout(spls_VarPlot, xaxis = list(title = xaxislabel,titlefont=list(size=fontsizes[1]),tickfont=list(size=fontsizes[2])),
                             yaxis = list(title = yaxislabel,titlefont=list(size=fontsizes[1]),tickfont=list(size=fontsizes[2])),
                             title = "sPLS Variable Plot",
                             titlefont = list(size=fontsizes[3]),
                             legend = list(font = list(size=fontsizes[4])))
      
      
      return(spls_VarPlot)
    }
    if(useannot==TRUE){
      #Generamos anotaciones
      
      annotations2 <- as.data.frame(row.names(vartable))
      colnames(annotations2)<-"IDENTIFIER"
      annotations <- left_join(annotations2,as.data.frame(annotations),by="IDENTIFIER") #con esto asignamos cada anotación a la variable correspondiente por rowname. Asi no hay fallo.
      annotations <- as.matrix(annotations) #Por alguna razon, sin esto, luego no se ordenan los factores ni los colores.
      rm(annotations2)
      
      #Sacamos los niveles del mapman, y las etiquetas
      mapmanlevels <- as.numeric(unique(annotations[order(as.numeric(annotations[,3])),3]))
      mapmannames <- annotations[order(as.numeric(annotations[,3])),]
      mapmannames <- as.vector(unique(mapmannames[,4]))
      
      spls_VarPlot <- plot_ly(as.data.frame(vartable), 
                              x = ~vartable[,1], 
                              y = ~vartable[,2],
                              text = paste("ID:", rownames(vartable),"</br> Description: ",annotations[,2],"</br> Mercator Bin: ",annotations[,4]), 
                              type="scatter",
                              mode = "markers", 
                              marker = list(size = 6, opacity=1),
                              color = factor(annotations[,4], levels=mapmannames),
                              colors = paleta_mapman[mapmanlevels],
                              symbol = factor(vartable[,3], labels=unique(vartable[,3])),
                              symbols = simbolos[1:length(unique(vartable[,3]))])
      
      spls_VarPlot <- layout(spls_VarPlot, xaxis = list(title = xaxislabel,titlefont=list(size=fontsizes[1]),tickfont=list(size=fontsizes[2])),
                             yaxis = list(title = yaxislabel,titlefont=list(size=fontsizes[1]),tickfont=list(size=fontsizes[2])),
                             title = "sPLS Variable Plot",
                             titlefont = list(size=fontsizes[3]),
                             legend = list(font = list(size=fontsizes[4])))
      
      return(spls_VarPlot)
    }
  }
  
  if(plottype=="CompoX-Y"){
    p1 <- spls_plot(spls_list,treatment, plottype = "Distance", fontsizes = fontsizes,variablespace = "x")
    p2 <- spls_plot(spls_list,treatment, plottype = "Distance", fontsizes = fontsizes,variablespace = "y")
    spls_CompositeplotPlot <- subplot(subplot(p1,p2),nrows=2,shareX = F,shareY = F,margin = 0.04)
    spls_CompositeplotPlot <- layout(spls_CompositeplotPlot, title = "SPLS Analysis",titlefont = list(size=fontsizes[3]))
    return(spls_CompositeplotPlot)
  }
  
  if(plottype=="Composite"){
    p1 <- spls_plot(spls_list,treatment, plottype = "Distance", fontsizes = fontsizes,variablespace = "x")
    p2 <- spls_plot(spls_list,treatment, plottype = "Distance", fontsizes = fontsizes,variablespace = "y")
    p3 <- spls_plot(spls_list,treatment, plottype = "Distance", fontsizes = fontsizes,variablespace = "xy")
    p4 <- spls_plot(spls_list,treatment, plottype = "Var", fontsizes = fontsizes)
    spls_CompositeplotPlot <- subplot(subplot(p1,p2),subplot(p3,p4),nrows=2,shareX = F,shareY = F,margin = 0.04)
    spls_CompositeplotPlot <- layout(spls_CompositeplotPlot, title = "SPLS Analysis",titlefont = list(size=fontsizes[3]))
    return(spls_CompositeplotPlot)
  }
  
  if(plottype=="Composite_Annotated"){
    p1 <- spls_plot(spls_list,treatment, plottype = "Distance", fontsizes = fontsizes,variablespace = "x")
    p2 <- spls_plot(spls_list,treatment, plottype = "Distance", fontsizes = fontsizes,variablespace = "y")
    p3 <- spls_plot(spls_list,treatment, plottype = "Distance", fontsizes = fontsizes,variablespace = "xy")
    p4 <- spls_plot(spls_list,treatment, plottype = "Var", fontsizes = fontsizes,useannot = T)
    spls_CompositeplotPlot <- subplot(subplot(p1,p2),subplot(p3,p4),nrows=2,shareX = F,shareY = F,margin = 0.04)
    spls_CompositeplotPlot <- layout(spls_CompositeplotPlot, title = "SPLS Analysis",titlefont = list(size=fontsizes[3]),showlegend=FALSE)
    return(spls_CompositeplotPlot)
  }
  
  if(plottype=="Network"){
    if(useannot==FALSE) {
      set.seed(5881)
      network(spls_object, comp=c(1,2),cutoff=cutoff)
      class(spls_object) <- c(class(spls_object),"networkplot")
      return(spls_object)
    }
    if(useannot==TRUE){
      set.seed(5881)
      network(spls_list[[5]], comp=c(1,2),cex.node.name = 0.5,cutoff=cutoff)
      class(spls_list[[5]]) <- c(class(spls_list[[5]]),"networkplotannot")
      return(spls_list[[5]])
    }
  }
}


