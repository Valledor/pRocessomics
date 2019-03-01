##
##  Diablo
##

## KeepX puede ser NULL - todas las variables; Auto - Se determina automaticamente la mejor opcion; Custom - se pasa una lista.
  
diablo_analysis <- function(listadatos, omiclevel=NULL,annotation=NULL, initialrow=1, initialcolumn=3, treatment1col=1, treatment2col=2, 
                            treatment=1, ncomponents=2, keepX=NULL,folds=3,nrepeat=5,autotune=FALSE,strength=0.1){
  options(warn=-1)
  require(mixOmics)
  require(plyr)
  require(doParallel)
  
  if(class(listadatos)!="POL") stop("A POL object is required. Please run pre-processing steps")
  datasetnames<-names(listadatos)
  
  if(is.null(omiclevel)==FALSE){
    sapply(omiclevel, function(x) if(x %in% datasetnames==FALSE) stop("Error: Please select a valid datasetname"))
    listadatos<-listadatos[omiclevel]
    datasetnames<-names(listadatos)
  }
  
  
  #Generamos un vector para hacer el split por tratamientos
  if(treatment%in%c(1,2,3)==FALSE) stop("Select a valid treatment")
  if(treatment==1) vectortratamientos <- as.factor(listadatos[[1]][initialrow:nrow(listadatos[[1]]),treatment1col])
  if(treatment==2) vectortratamientos <- as.factor(listadatos[[1]][initialrow:nrow(listadatos[[1]]),treatment2col])
  if(treatment==3) vectortratamientos <- as.factor(paste(listadatos[[1]][initialrow:nrow(listadatos[[1]]),treatment1col],listadatos[[1]][initialrow:nrow(listadatos[[1]]),treatment2col], sep="&"))
  
  #aprovechamos para hacer la tabla de tratamientos. si no se hace esto los plots pierden versatilidad y o petan
  if(is.null(treatment2col)) treatment2col<-treatment1col
  if(treatment1col==treatment2col) tratamientos <-cbind(listadatos[[1]][initialrow:nrow(listadatos[[1]]),treatment1col],listadatos[[1]][initialrow:nrow(listadatos[[1]]),treatment1col])
  if(treatment1col!=treatment2col) tratamientos <-cbind(listadatos[[1]][initialrow:nrow(listadatos[[1]]),treatment1col],listadatos[[1]][initialrow:nrow(listadatos[[1]]),treatment2col])
  tratamientos <- cbind(tratamientos, paste(listadatos[[1]][initialrow:nrow(listadatos[[1]]),treatment1col],listadatos[[1]][initialrow:nrow(listadatos[[1]]),treatment2col], sep="&"))
  
  #Extraemos los valores significativos de las tablas
  listadatos2 <- lapply(listadatos, function(x) x[initialrow:nrow(x),initialcolumn:ncol(x)])
  names(listadatos2)<-names(listadatos)
  
  #Miramos que las dimensiones sean correctas. X ha de ser igual para todas las matrices
  dimlist<-unlist(lapply(lapply(listadatos2,dim),function(x) x[1]))
  if(length(unique(dimlist))>1) stop("Cases in all datasets should be the same. Please re-check nrows")
  
  #Diseñamos la matriz
  matrixdesign <- matrix(strength, ncol=length(listadatos2), nrow =length(listadatos2), dimnames = list(datasetnames,datasetnames))
  diag(matrixdesign) <-0
  
  if(!is.null(annotation)){ #Creamos una tabla de anotaciones, a nuestra manera, para plotear luego. Chapuza, pero funciona.
    annotation_table <- annotation_mapman(listadatos2,annotation)
    #ahora hacemos lo mismo pero en formato lista para poder poner nombre a las redes exportadas
    annotation_list <- lapply(listadatos2, function(x) annotation_mapman(as.data.frame(x),annotation))
    names(annotation_list)<-names(listadatos)
  }
  
  if(autotune==TRUE){
    options(warn=-1)
    #Dos componentes está bien. Ahora evaluamos el número de variables a mantener
    pruebaX <-lapply(listadatos2,series)
    cat("\nTuning DIABLO parameters. It may take a while")
    library(parallel)
    cores<-detectCores()
    set.seed(5881)
    tune.diablo <- tune.block.splsda(X=listadatos2,Y=vectortratamientos,ncomp=ncomponents,test.keepX = pruebaX, design = matrixdesign, validation = "Mfold", folds=folds, nrepeat=nrepeat, dist="centroids.dist",cpus=cores)
    
    #Creamos la lista de variables a mantener
    listakeepX <- tune.diablo$choice.keepX
    keepX <- listakeepX
    cat("Best errors in components achieved for:\n")
    print(keepX)
    readline(prompt="Press [enter] to continue")
  }
  
  #if(class(keepX)!="list") stop("Variables for KeepX should be provided as list")
  #if(Reduce("&",datasetnames==names(keepX))==FALSE) stop("List elements for KeepX should have the same name than elements in data list")
  #if(Reduce("&",sapply(listakeepX,class)=="numeric")==FALSE) stop("Elements in KeepX List must be numeric")
  #if(Reduce("&",sapply(listakeepX,length)==rep(ncomponents,length(datasetnames))==FALSE)) stop("Elements in the list must have the same number of elements than components")
  listakeepX <-keepX
  
  # Hacemos el análisis 
    cat("\nCALCULATING DIABLO\nIt may take a while...\n")
    diablo_analysis <- block.splsda(X=listadatos2,Y=vectortratamientos,ncomp=ncomponents,design = matrixdesign,keepX = listakeepX)
    
    if(is.null(annotation)){
      resultado <- list(diablo_analysis,as.data.frame(tratamientos),names(listadatos))
      names(resultado) <- c("diablo","treatments","datasetnames")
      class(resultado) <- "diabloanalysis"
    } else {
      #Esto es una chapuza para que las redes puedan salir en pantalla con nombres de variables
      #Este bucle de abajo da problemas si hay nombres similares en dos niveles omicos distintos
      #for(i in 1:length(listadatos2)){
      #  colnames(listadatos2[[i]])<-make.names(annotation_list[[i]][,2],unique=TRUE)
      #}
      #Alternativa
      vectorncolsets <- c()
      nombresvariables <- c()
      for(i in 1:length(listadatos2)){
        vectorncolsets<-c(vectorncolsets, ncol(listadatos2[[i]]))
        nombresvariables <- c(nombresvariables, annotation_list[[i]][,2])
      }
      nombresvariables <- make.names(nombresvariables,unique=TRUE)
      for(j in 1:length(vectorncolsets)){
        if(j==1) colnames(listadatos2[[j]])<-nombresvariables[1:vectorncolsets[j]]
        if(j>1) colnames(listadatos2[[j]])<-nombresvariables[sum(vectorncolsets[1:(j-1)]+1):sum(vectorncolsets[1:j])]
      }
      
      diablo_network<-block.splsda(X=listadatos2,Y=vectortratamientos,ncomp=ncomponents,design = matrixdesign,keepX = listakeepX)
      resultado <- list(diablo_analysis,as.data.frame(tratamientos),names(listadatos),annotation_table,diablo_network)
      names(resultado) <- c("diablo","treatments","datasetnames","annotations","diablo_network")
      class(resultado) <- "diabloanalysis"
      return(resultado)
    }
}

diablo_plot <- function(diablo_list,treatment=1, plottype="Distance",variablespace="xy",useannot=FALSE,fontsizes=c(14,10,16,12),cutoff=0.5){ 
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
  if(class(diablo_list)!= "diabloanalysis")
    stop("diablo analysis object is expected, please run diablo_analysis first")
  if(plottype %in% c("Composite","Distance", "Var","Network","All","Composite_Annotated","cim","Circos")==FALSE)
    stop("Please select a valid plot type")
  if(variablespace%in%c("x","y","xy")==FALSE) stop("Select an adequate variable representation space. x, y, or xy ")
  if(useannot==T&&is.null(diablo_list$annotations)){
    cat(paste("\nAnnotation matrix not provided. Considering useannot=F\n", "Please check/repeat da_analysis() including an annotation matrix\n","Read help or follow the tutorial for more information.\n",sep=""))
    useannot=FALSE
  }
  #...
  if(dev.cur()>1) { dev.off()}
  ##
  ## Definimos los objetos necesarios para dibujar
  ##
  
  diablo_object <- diablo_list$diablo
  treatmentstable <- diablo_list$treatments
  annotations <- diablo_list$annotations
  datasetnames <- diablo_list$datasetnames
  treatmentnames <- colnames(treatmentstable)
  
  #Plot sencillo, Distancias
  if(plottype=="Distance"){
    a<-list()
    for(i in 1:length(diablo_object$variates)){
      a[[i]]<-distanceplotdiablo(diablo_object$variates[[i]],diablo_object$explained_variance[[i]],treatment,treatmentstable,names(diablo_object$explained_variance)[i])  
    }
    diablo_distance_plot <- subplot(a,nrows = 2)
    diablo_distance_plot <- layout(diablo_distance_plot,
                                   title = paste("Distance Plots:",paste(names(diablo_object$explained_variance),collapse=", "),sep=" "),
                                   titlefont = list(size=fontsizes[3]),
                                   legend = list(font = list(size=fontsizes[4])))
    return(diablo_distance_plot)
  }
  
  if(plottype=="Var"){
    vartable <-invisible(plotVar(diablo_object,plot=F)) #esto es una chapuza. habria que diseñar funcion para sacar esto. de momento uso mixomics
    xaxislabel <- paste("Component 1")
    yaxislabel <- paste("Component 2")
    
    if(useannot==FALSE){
      diablo_VarPlot <- plot_ly(as.data.frame(vartable), 
                                x = ~vartable[,1], 
                                y = ~vartable[,2],
                                text = paste("ID:", rownames(vartable)), 
                                type="scatter",
                                mode = "markers", 
                                marker = list(size = 6, opacity=1), 
                                colors =paleta_niveles[1:length(datasetnames)],
                                symbol = factor(vartable[,3], labels=unique(vartable[,3])),
                                symbols = simbolos[1:length(unique(vartable[,3]))])
      
      diabloVarPlot <- layout(diablo_VarPlot, xaxis = list(title = xaxislabel,titlefont=list(size=fontsizes[1]),tickfont=list(size=fontsizes[2])),
                              yaxis = list(title = yaxislabel,titlefont=list(size=fontsizes[1]),tickfont=list(size=fontsizes[2])),
                              title = "Diablo Variable Plot",
                              titlefont = list(size=fontsizes[3]),
                              legend = list(font = list(size=fontsizes[4])))
      
      
      return(diablo_VarPlot)
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
      
      diablo_VarPlot <- plot_ly(as.data.frame(vartable), 
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
      
      diablo_VarPlot <- layout(diablo_VarPlot, xaxis = list(title = xaxislabel,titlefont=list(size=fontsizes[1]),tickfont=list(size=fontsizes[2])),
                               yaxis = list(title = yaxislabel,titlefont=list(size=fontsizes[1]),tickfont=list(size=fontsizes[2])),
                               title = "Diablo Variable Plot",
                               titlefont = list(size=fontsizes[3]),
                               legend = list(font = list(size=fontsizes[4])))
      
      return(diablo_VarPlot)
    }
  }
  
  if(plottype=="cim"){
    if(useannot==T){
    diablo_object <-diablo_list[[5]]
    cimDiablo(diablo_object,legend.position = 0,size.legend=1,margins = c(8,8))
    class(diablo_object) <- c(class(diablo_object),"cimdiablo")
    return(diablo_object)
    } else {
      cimDiablo(diablo_object,legend.position = 0,size.legend=1,margins = c(8,8))
      class(diablo_object) <- c(class(diablo_object),"cimdiablo")
      return(diablo_object)
    }
    
    
  }
  
  if(plottype=="Network"){
    if(useannot==FALSE) {
      set.seed(5881)
      network(diablo_object, cutoff=cutoff)
      class(diablo_object) <- c(class(diablo_object),"networkplot")
      return(diablo_object)
    }
    if(useannot==TRUE){
      set.seed(5881)
      network(diablo_list[[5]],cex.node.name = 0.5,cutoff=cutoff)
      class(diablo_list[[5]]) <- c(class(diablo_list[[5]]),"networkplotannot")
      return(diablo_list[[5]])
    }
  }
  
  if(plottype=="Circos"){
    if(useannot==FALSE) {
      set.seed(5881)
      circosPlot(diablo_object, cutoff=cutoff)
      class(diablo_object) <- c(class(diablo_object),"circosplot")
      return(diablo_object)
    }
    if(useannot==TRUE){
      set.seed(5881)
      circosPlot(diablo_list[[5]],cutoff=cutoff)
      class(diablo_list[[5]]) <- c(class(diablo_list[[5]]),"circosplot")
      return(diablo_list[[5]])
    }
  }
}

distanceplotdiablo<- function(variates,varianceexplained,treatment, treatmentstable,levelname,fontsizes=c(14,10,16,12)){
  explainedvariance <-as.vector(varianceexplained)
  xaxislabel <- paste("Component 1-",round(explainedvariance[1]*100,2),"% of variance", sep=" ")
  yaxislabel <- paste("Component 2-",round(explainedvariance[2]*100,2),"% of variance", sep=" ")
  
  if(treatment!=4){  #Asi seleccionamos todas las visualizaciones que solo tengan un factor a mostrar.
    if(treatment==1) tratamientocolores <-as.vector(treatmentstable[,1])
    if(treatment==2) tratamientocolores <-as.vector(treatmentstable[,2])
    if(treatment==3) tratamientocolores <-as.vector(treatmentstable[,3])
    
    da_distance_plot <- plot_ly(as.data.frame(variates), 
                                x = ~variates[,1], 
                                y = ~variates[,2],
                                hoverinfo ="text",
                                text = paste("Sample:", rownames(variates), sep=" "),
                                showlegend = TRUE,
                                type="scatter",
                                mode = "markers", 
                                color = factor(tratamientocolores,labels = c(unique(tratamientocolores))),
                                colors = paleta_tratamientos[1:length(unique(tratamientocolores))],
                                #symbol = factor(rep(c(1:nrow(mcoin_object$mcoa$lambda)),each=nrow(mcoin_object$mcoa$SynVar)), labels=datasetnames),
                                #symbols = simbolos[1:length(unique(tratamientocolores))], 
                                marker = list(size = 11)) %>%
      layout(xaxis = list(title = xaxislabel,titlefont=list(size=fontsizes[1]),tickfont=list(size=fontsizes[2])),
             yaxis = list(title = yaxislabel,titlefont=list(size=fontsizes[1]),tickfont=list(size=fontsizes[2])),
             title = levelname,
             titlefont = list(size=fontsizes[3]),
             legend = list(font = list(size=fontsizes[4])))
    
    return(da_distance_plot)
  }
  
  if(treatment==4){  #Asi seleccionamos las visualizaciones con 2 factores a mostrar
    
    da_distance_plot <- plot_ly(as.data.frame(variates), 
                                x = ~variates[,1], 
                                y = ~variates[,2],
                                hoverinfo ="text",
                                text = paste("Sample:", rownames(variates), sep=" "),
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
             title = levelname,
             titlefont = list(size=fontsizes[3]),
             legend = list(font = list(size=fontsizes[4])))
    return(da_distance_plot)
    
  }
}
