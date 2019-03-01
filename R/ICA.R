ica_analysis <- function(listadatos, annotation=NULL, initialrow=1, initialcolumn=3,  treatment1col=1, treatment2col=2, ncomp=3,omiclevel="all"){
  options(warn=-1)
  require("ica")
  if(class(listadatos)!="POL") stop("A POL object is required. Please run pre-processing steps")
  datasetnames<-names(listadatos)
  
  #Extraemos los valores significativos de las tablas
  listadatos2 <- lapply(listadatos, function(x) x[initialrow:nrow(x),initialcolumn:ncol(x)])
  if(is.null(treatment2col)) treatment2col=treatment1col
  
  # Seleccionamos los datos
  if(is.null(omiclevel)|omiclevel=="all") { #esto no va!
    variablenames <- lapply(listadatos2, function(x) colnames(x))
    variablenames <-unlist(variablenames)
    listaintermedia <-do.call(cbind,listadatos2)  #### Esto es inconsistente. Cambia nombres de columnas, a veces, nombretabla.nombrevar, otras veces solo nombrevar
    colnames(listaintermedia) <- variablenames #### para arreglar lode arriba entre tanto.
    listaintermedia<-RemoveEmptyColumns(listaintermedia,1,1)#Esto es para evitar cosas raras.
  }  
  if(omiclevel!="all"&!is.null(omiclevel)) {
    sapply(omiclevel, function(x) if(x %in% datasetnames==FALSE) stop("Error: Please select a valid datasetname"))
    
    if(length(omiclevel)>1){  listaintermedia <-do.call(cbind,listadatos2[omiclevel])
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
  
  ##
  ## Apano para crear una tabla de identificaciones a usar posteriormente. Solo MapMan.
  ## es un poco chapuzas pero bueno, funciona
  ##
  #
  # Importante: chequear que la tabla de identificaciones no tenga duplicados, sino es un lio porque el dplyr crea filas extra y no cuadra 
  # para los dibujos. Pensar en algo para que si hay duplicados de error y no siga.
  
  if(!is.null(annotation)){ #Creamos una tabla de anotaciones, a nuestra manera, para plotear luego. Chapuza, pero funciona.
    library(dplyr)
    columnnames <- c("IDENTIFIER","DESCRIPTION","MAPMAN.BINCODE","MAPMAN.BINNAME")
    annotmatrix <-cbind(variablenames,"") # esto es una caca, pero sino el dplyr no me funciona
    colnames(annotmatrix) <- c("IDENTIFIER","")
    annotmatrix <- left_join(as.data.frame(annotmatrix),annotation,by="IDENTIFIER")
    annotmatrix <- annotmatrix[,-2]
    annotmatrix <- sapply(annotmatrix,as.vector)
    annotmatrix[annotmatrix[,2]=="",2]<-"" #Dejamos en blanco todo aquello que no tenga anotaciones
    annotmatrix[annotmatrix[,3]=="",3]<-35.2 #Definimos como desconocidos los valores perdidos
    annotmatrix[annotmatrix[,4]=="",4]<-"not assigned" # lo mismo.
    annotmatrix[,3] <- as.numeric(do.call(rbind,lapply(strsplit(as.vector(annotmatrix[,3]),"\\."), function(x) x[[1]]))) #me quedo con el primer bin
    annotmatrix[,4] <- as.character(do.call(rbind,lapply(strsplit(as.vector(annotmatrix[,4]),"\\."), function(x) x[[1]]))) #me quedo con el primer descriptor
    colnames(annotmatrix) <- columnnames
    annotmatrix[,2][is.na(annotmatrix[,2])] <-""
    annotmatrix[,3][is.na(annotmatrix[,3])] <-"35"
    annotmatrix[,4][is.na(annotmatrix[,4])] <-"not assigned"
    annotation<-annotmatrix
  }
  
  # Definimos tratamientos
  if(treatment1col==treatment2col){
    treatments <- as.data.frame(as.vector(listadatos[[1]][initialrow:nrow(listadatos[[1]]),treatment1col]))
    colnames(treatments) <-colnames(listadatos[[1]][,treatment1col])
    #treatments <- cbind(treatments,treatments)
    
  } else {
    treatments<-listadatos[[1]][initialrow:nrow(listadatos[[1]]),c(treatment1col:treatment2col)]
  }
  
  # Hacemos el análisis 
  cat("\nCALCULATING ICA\nIt may take a while...\n")
  ica_analysis <- icafast(listaintermedia,nc=ncomp,alg="def",center=T,fun="logcosh")
  colnames(ica_analysis$W)<-colnames(listaintermedia)
  
  if(is.null(annotation)){
    resultado <- list(ica_analysis,treatments,names(listadatos))
    names(resultado) <- c("ica","treatments","datasetnames")
    class(resultado) <- "icaanalysis"
    return(resultado)
  } else {
    resultado <- list(ica_analysis,treatments,names(listadatos),annotation)
    names(resultado) <- c("ica","treatments","datasetnames","annotations")
    class(resultado) <- "icaanalysis"
    return(resultado)
  }
}

ica_plot <- function(ica_list,treatment=1, compX=1, compY=2, compZ=NULL, plottype="Composite",useannot=FALSE){ 
  library("plotly") 
  ## Se ha puesto un tipo de plot por defecto porque si no, no sale nada, pero no devuelve error
  ## Se dara la opcion de escoger un tratamiento, que sera factor 1[1], factor 2 [2](si lo definio el usuario), interaccion[3], o combinar ambos[4].
  ## Falta corregir algunas funciones graficas, y temas de maquetado/leyendas, pero es cosa menor y no afecta a la funcionalidad
  ## Es guay, se puede llamar a una funcion de forma recursiva desde la misma funcion... (esto no lo sabia)
  
  ## Hay un bloque de vectores de "imagen" con los sets de colores y de simbolos. Como habra que cambiarlos, porque no seran
  ## todos del agrado de la mayoria, pues mejor que se definan fuera de todo.
  
  ## Plottype: Aqui estan, o estaran, todos los graficos que pueden ser utiles. Hay unos cuantos, la gente lega no los usara, 
  ## pero es interesante mostrarlos.
  #             "Composite" <- Collage con los 2 plots mas representativos del analisis
  #             "Score" <- representa los individuos (1 punto por muestra), considerando todos los niveles. Unico que permitira 3D
  #             "Scree" <- Scree plot similar al del PCA
  #             "Var" <- Muestrea las correlaciones de variables, permite escoger dataset y anotacion
  #             "All" <- Devuelve una lista con todos los graficos anteriores.
  
  
  ## ********Esto hay que mejorarlo un poco controlando todos los niveles del objeto y el tipo de grafico*******
  ## Validamos que el objeto de verdad sea del tipo que tiene que ser:
  ##
  if(class(ica_list)!= "icaanalysis")
    stop("ICA analysis object is expected, please run ica_analysis first")
  if(!inherits(ica_list[[1]], "list"))
    stop("ICA analysis object seems to be wrong. Please re-run ica_analysis")
  if(plottype %in% c("Composite","Score","Scree","Var","All")==FALSE)
    stop("Please select a valid plot type")
  if(useannot==T&&is.null(ica_list$annotations)){
    cat(paste("\nAnnotation matrix not provided. Considering useannot=F\n", "Please check/repeat ica_analysis() including an annotation matrix\n","Read help or follow the tutorial for more information.\n",sep=""))
    useannot=FALSE
  }
  #...
  
  ##
  ## Definimos los objetos necesarios para dibujar
  ##
  
  ica_object <- ica_list$ica
  treatmentstable <- ica_list$treatments
  annotations <- ica_list$annotations
  datasetnames <- ica_list$datasetnames
  treatmentnames <- colnames(treatmentstable)
  
  #Plot sencillo, Scores
  if(plottype=="Score"){
    if(compX>length(ica_object$vafs)) stop(paste("\nX-Axis component do not exist. Only ",length(ica_object$vafs)," components were defined in ICA analysis",sep=""))
    if(compY>length(ica_object$vafs)) stop(paste("\nY-Axis component do not exist. Only ",length(ica_object$vafs)," components were defined in ICA analysis",sep=""))
    
    if(treatment!=3){  #Asi seleccionamos todas las visualizaciones que solo tengan un factor a mostrar.
      if(treatment==1) tratamientocolores <-as.vector(treatmentstable[,1])
      if(treatment==2) tratamientocolores <-as.vector(treatmentstable[,2])
      if(treatment==4) tratamientocolores <- as.vector(paste(treatmentstable[,1],treatmentstable[,2],sep="&"))
      if(is.null(compZ)){ # Plot 2D, 1 factor.
        
        ica_score_plot <- plot_ly(as.data.frame(ica_object$S), 
                                  x = ~ica_object$S[,compX], 
                                  y = ~ica_object$S[,compY],
                                  hoverinfo ="text",
                                  text = paste("Sample:", rownames(ica_object$S), sep=" "),
                                  showlegend = TRUE,
                                  type="scatter",
                                  mode = "markers", 
                                  color = factor(tratamientocolores,labels = c(unique(tratamientocolores))),
                                  colors = paleta_tratamientos[1:length(unique(tratamientocolores))],
                                  #symbol = factor(rep(c(1:nrow(mcoin_object$mcoa$lambda)),each=nrow(mcoin_object$mcoa$SynVar)), labels=datasetnames),
                                  symbols = simbolos[1:length(unique(tratamientocolores))], 
                                  marker = list(size = 11)) %>%
          layout(xaxis = list(title = paste("Component",compX,"-",round(ica_object$vafs*100,2)[compX],"% of variance", sep=" "),titlefont=list(size=fontsizes[1]),tickfont=list(size=fontsizes[2])),
                 yaxis = list(title = paste("Component",compY,"-",round(ica_object$vafs*100,2)[compY],"% of variance",sep=" "),titlefont=list(size=fontsizes[1]),tickfont=list(size=fontsizes[2])),
                 title = "ICA ScorePlot",
                 titlefont = list(size=fontsizes[3]),
                 legend = list(font = list(size=fontsizes[4])))
        
        return(ica_score_plot)
      }
      
      if(!is.null(compZ)){ # Plot 3D, 1 factor.
        if(compZ>length(ica_object$vafs)) stop(paste("\nZ-Axis component do not exist. Only ",length(ica_object$vafs)," components were defined in ICA analysis",sep=""))
        ica_score_plot <- plot_ly(as.data.frame(ica_object$S), 
                                  x = ~ica_object$S[,compX], 
                                  y = ~ica_object$S[,compY],
                                  z = ~ica_object$S[,compZ],
                                  hoverinfo ="text",
                                  text = paste("Sample:", rownames(ica_object$S), sep=" "),
                                  showlegend = TRUE,
                                  type="scatter3d",
                                  mode = "markers", 
                                  color = factor(tratamientocolores,labels = c(unique(tratamientocolores))),
                                  colors = paleta_tratamientos[1:length(unique(tratamientocolores))],
                                  
                                  symbols = simbolos[1:length(unique(tratamientocolores))], 
                                  marker = list(size = 11)) %>%
          
          layout(scene=list(xaxis = list(title = paste("Component",compX,"-",round(ica_object$vafs*100,2)[compX],"% of variance", sep=" "),titlefont=list(size=fontsizes[1]),tickfont=list(size=fontsizes[2])),
                 yaxis = list(title = paste("Component",compY,"-",round(ica_object$vafs*100,2)[compY],"% of variance",sep=" "),titlefont=list(size=fontsizes[1]),tickfont=list(size=fontsizes[2])),
                 zaxis = list(title = paste("Component",compZ,"-",round(ica_object$vafs*100,2)[compZ],"% of variance",sep=" "),titlefont=list(size=fontsizes[1]),tickfont=list(size=fontsizes[2]))),
                 title = "ICA ScorePlot",
                 titlefont = list(size=fontsizes[3]),
                 legend = list(font = list(size=fontsizes[4])))
        
        return(ica_score_plot)
      }
      
    }
    if(treatment==3){  #Asi seleccionamos las visualizaciones con 2 factores a mostrar
      
      if(is.null(compZ)){ # Plot 2D, 2 factores.
        ica_score_plot <- plot_ly(as.data.frame(ica_object$S), 
                                  x = ~ica_object$S[,compX], 
                                  y = ~ica_object$S[,compY],
                                  hoverinfo ="text",
                                  text = paste("Sample:", rownames(ica_object$S), sep=" "),
                                  showlegend = TRUE,
                                  type="scatter",
                                  mode = "markers", 
                                  color = factor(treatmentstable[,1], labels = unique(treatmentstable[,1])),
                                  colors = paleta_tratamientos[1:length(unique(treatmentstable[,1]))],
                                  symbol = factor(treatmentstable[,2], labels=unique(treatmentstable[,2])),
                                  symbols = simbolos[1:length(unique(treatmentstable[,2]))], 
                                  marker = list(size = 11)) %>%
          layout(xaxis = list(title = paste("Component",compX,"-",round(ica_object$vafs*100,2)[compX],"% of variance", sep=" "),titlefont=list(size=fontsizes[1]),tickfont=list(size=fontsizes[2])),
                 yaxis = list(title = paste("Component",compY,"-",round(ica_object$vafs*100,2)[compY],"% of variance",sep=" "),titlefont=list(size=fontsizes[1]),tickfont=list(size=fontsizes[2])),
                 title = "ICA ScorePlot",
                 titlefont = list(size=fontsizes[3]),
                 legend = list(font = list(size=fontsizes[4])))
        return(ica_score_plot)
      }
      
      if(!is.null(compZ)){ # Plot 3D, 2 factores.
        if(compZ>length(ica_object$vafs)) stop(paste("\nZ-Axis component do not exist. Only ",length(ica_object$vafs)," components were defined in ICA analysis",sep=""))
        
        ica_score_plot <- plot_ly(as.data.frame(ica_object$S), 
                                  x = ~ica_object$S[,compX], 
                                  y = ~ica_object$S[,compY],
                                  z = ~ica_object$S[,compZ],
                                  hoverinfo ="text",
                                  text = paste("Sample:", rownames(ica_object$S), sep=" "),
                                  showlegend = TRUE,
                                  type="scatter3d",
                                  mode = "markers", 
                                  color = factor(treatmentstable[,1], labels = unique(treatmentstable[,1])),
                                  colors = paleta_tratamientos[1:length(unique(treatmentstable[,1]))],
                                  symbol = factor(treatmentstable[,2], labels=unique(treatmentstable[,2])),
                                  symbols = simbolos[1:length(unique(treatmentstable[,2]))], 
                                  marker = list(size = 11)) %>%
          layout(scene=list(xaxis = list(title = paste("Component",compX,"-",round(ica_object$vafs*100,2)[compX],"% of variance", sep=" "),titlefont=list(size=fontsizes[1]),tickfont=list(size=fontsizes[2])),
                            yaxis = list(title = paste("Component",compY,"-",round(ica_object$vafs*100,2)[compY],"% of variance",sep=" "),titlefont=list(size=fontsizes[1]),tickfont=list(size=fontsizes[2])),
                            zaxis = list(title = paste("Component",compZ,"-",round(ica_object$vafs*100,2)[compZ],"% of variance",sep=" "),titlefont=list(size=fontsizes[1]),tickfont=list(size=fontsizes[2]))),
                 title = "ICA ScorePlot",
                 titlefont = list(size=fontsizes[3]),
                 legend = list(font = list(size=fontsizes[4])))
        
        return(ica_score_plot)
        
      }
    }
  }
  
  if(plottype=="Scree"){
    screeplotdata <- data.frame(ICs=paste("IC",c(1:length(ica_object$vafs)),sep=""), ExplVar=round(ica_object$vafs*100,2), Cumulat=round(cumsum(ica_object$vafs*100),2),stringsAsFactors = FALSE) #Creamos tabla de datos con % varianza y % var acumulado
    screeplotdata$ICs <- factor(screeplotdata$ICs, levels = unique(screeplotdata$ICs)[order(screeplotdata$ExplVar, decreasing = TRUE)]) # Ordenamos los ICs segun valor
    
    ica_screeplot <- plot_ly(screeplotdata) %>%
      add_trace(x=~ICs,y=~ExplVar,type="bar",name="% Explained Variance",
                marker = list(color = '#C9EFF9')) %>%
      add_trace(x=~ICs, y=~Cumulat,type="scatter",mode="lines+markers",yaxis = "y2",name="Cumulative Proportion",
                line = list(color = '#45171D')) %>%
    layout(title = "ICA Screeplot",
           titlefont = list(size=fontsizes[3]),
           xaxis = list(title = "",titlefont=list(size=fontsizes[1]),tickfont=list(size=fontsizes[2])),
           yaxis = list(side = 'left', title = "% of Explained Variance", showgrid = F, zeroline = TRUE,titlefont=list(size=fontsizes[1]),tickfont=list(size=fontsizes[2])),
           yaxis2 = list(side = 'right', overlaying = "y", title = "Cumulative", showgrid = FALSE, zeroline = FALSE,titlefont=list(size=fontsizes[1]),tickfont=list(size=fontsizes[2])),
           legend = list(x = 0.7, y = 0.1),font = list(size=fontsizes[4]))
      return(ica_screeplot) 
    
  }
  
  if(plottype=="Var"){
    #Generamos anotaciones
    if(!is.null(compZ)) cat("\nVarplot can only plot 2 components. Only plotting X and Y.")
    annotations2 <- as.data.frame(colnames(ica_object$W))
    colnames(annotations2)<-"IDENTIFIER"
    annotations <- left_join(annotations2,as.data.frame(annotations),by="IDENTIFIER") #con esto asignamos cada anotación a la variable correspondiente por rowname. Asi no hay fallo.
    annotations <- as.matrix(annotations) #Por alguna razon, sin esto, luego no se ordenan los factores ni los colores.
    rm(annotations2)
    
    if(useannot==FALSE){
      ica_VarPlot <- plot_ly(as.data.frame(ica_object$W), 
                             x = ~ica_object$W[compX,], 
                             y = ~ica_object$W[compY,],
                             text = paste("ID:", colnames(ica_object$W)), 
                             type="scatter",
                             mode = "markers", 
                             marker = list(size = 6, opacity=1), 
                             colors =paleta_niveles[1:length(datasetnames)])

      pca_VarPlot <- layout(ica_VarPlot,
                            xaxis = list(title = paste("Component",compX,"-",explainedvariance[compX],"% of variance", sep=" "),titlefont=list(size=fontsizes[1]),tickfont=list(size=fontsizes[2])),
                            yaxis = list(title = paste("Component",compY,"-",explainedvariance[compY],"% of variance",sep=" "),titlefont=list(size=fontsizes[1]),tickfont=list(size=fontsizes[2])),
                            title = "ICA Estimated Unmixing Signals",
                            titlefont = list(size=fontsizes[3]),
                            legend = list(font = list(size=fontsizes[4]))) 
      
      
      return(ica_VarPlot)
    }
    
    if(useannot==TRUE){
      #Sacamos los niveles del mapman, y las etiquetas
      mapmanlevels <- as.numeric(unique(annotations[order(as.numeric(annotations[,3])),3]))
      mapmannames <- annotations[order(as.numeric(annotations[,3])),]
      mapmannames <- as.vector(unique(mapmannames[,4]))
      
      ica_VarPlot <- plot_ly(as.data.frame(ica_object$W), 
                             x = ~ica_object$W[compX,], 
                             y = ~ica_object$W[compY,],
                             text = paste("ID:", colnames(ica_object$W),"</br> Description: ",annotations[,2],"</br> Mercator Bin: ",annotations[,4]), 
                             type="scatter",
                             mode = "markers", 
                             marker = list(size = 6, opacity=1),
                             color = factor(annotations[,4], levels=mapmannames),
                             colors = paleta_mapman[mapmanlevels])
      
      pca_VarPlot <- layout(ica_VarPlot,
                            xaxis = list(title = paste("Component",compX,"-",explainedvariance[compX],"% of variance", sep=" "),titlefont=list(size=fontsizes[1]),tickfont=list(size=fontsizes[2])),
                            yaxis = list(title = paste("Component",compY,"-",explainedvariance[compY],"% of variance",sep=" "),titlefont=list(size=fontsizes[1]),tickfont=list(size=fontsizes[2])),
                            title = "ICA Estimated Unmixing Signals",
                            titlefont = list(size=fontsizes[3]),
                            legend = list(font = list(size=fontsizes[4]))) 
      
      return(ica_VarPlot)
    }
  }
  
  if(plottype=="Composite"){
    if(!is.null(compZ)) stop("Composite plot can only draw 2D plots. compZ must be null") 
    #Representamos los eigenvalues, y luego como sube la proporcion acumulada.
    p1 <- ica_plot(ica_list,treatment, compX, compY, plottype = "Score")
    p2 <- ica_plot(ica_list,treatment, compX, compY, plottype = "Var", useannot=useannot)
    ica_composite <- subplot(p1,p2,nrows=2,shareX = F,shareY = F,margin = 0.04)
    ica_composite <- layout(ica_composite, title = "ICA Analysis")
    return(ica_composite)
  }
  
  if(plottype=="All"){
    #Aqui habra que crear una llamada recursiva a todos los graficos posibles y devolver una lista.
    #tambien recursivo pa los distintos niveles.
    if(!is.null(compZ)) stop("Please select only two dimensions for plotting all figures")
    cat("EXPORTING ICA ANALYSIS PLOTS\nA single core will be used. It may take a while...\n ")
    cat("Generating and saving Score plot. Step 1/4\n")
    export_plot(ica_plot(ica_list,treatment=treatment, compX=compX, compY=compY, plottype="Score",useannot=useannot),"ica_score.pdf")
    cat("Generating and saving Scree plot. Step 2/4\n")
    export_plot(ica_plot(ica_list,treatment=treatment, compX=compX, compY=compY, plottype="Scree",useannot=useannot),"ica_scree.pdf")
    cat("Generating and saving Var plot. Step 3/4\n")
    export_plot(ica_plot(ica_list,treatment=treatment, compX=compX, compY=compY, plottype="Var",useannot=useannot),"ica_var.pdf")
    cat("Generating and saving Composite plot. Step 4/4\n")
    export_plot(ica_plot(ica_list,treatment=treatment, compX=compX, compY=compY, plottype="Composite",useannot=useannot),"ica_composite.pdf")
  }
}