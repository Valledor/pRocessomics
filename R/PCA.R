pca_analysis <- function(listadatos, annotation=NULL, initialrow=1, initialcolumn=3, treatment1col=1, treatment2col=2, treatment=1,omiclevel="all"){
  options(warn=-1)
  if(class(listadatos)!="POL") stop("A POL object is required. Please run pre-processing steps")
  datasetnames<-names(listadatos)

  #Generamos un vector para hacer el split por tratamientos
  if(treatment%in%c(1,2,3)==FALSE) stop("Select a valid treatment")
  if(treatment==1) vectortratamientos <- as.factor(listadatos[[1]][initialrow:nrow(listadatos[[1]]),treatment1col])
  if(treatment==2) vectortratamientos <- as.factor(listadatos[[1]][initialrow:nrow(listadatos[[1]]),treatment2col])
  if(treatment==3) vectortratamientos <- as.factor(paste(listadatos[[1]][initialrow:nrow(listadatos[[1]]),treatment1col],listadatos[[1]][initialrow:nrow(listadatos[[1]]),treatment2col], sep="&"))

  #aprovechamos para hacer la tabla de tratamientos. si no se hace esto los plots pierden versatilidad y o petan
  if(is.null(treatment2col)) treatment2col<-treatment1col
  if(treatment1col==treatment2col) tratamientos <-cbind(as.character(listadatos[[1]][initialrow:nrow(listadatos[[1]]),treatment1col]),as.character(listadatos[[1]][initialrow:nrow(listadatos[[1]]),treatment1col]))
  if(treatment1col!=treatment2col) tratamientos <-cbind(as.character(listadatos[[1]][initialrow:nrow(listadatos[[1]]),treatment1col]),as.character(listadatos[[1]][initialrow:nrow(listadatos[[1]]),treatment2col]))
  tratamientos <- cbind(tratamientos, paste(listadatos[[1]][initialrow:nrow(listadatos[[1]]),treatment1col],listadatos[[1]][initialrow:nrow(listadatos[[1]]),treatment2col], sep="&"))

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

  if(!is.null(annotation)){
      annotation_table <- annotation_mapman(listaintermedia,annotation)
  }

  # Hacemos el análisis
  cat("\nCALCULATING PCA\nIt may take a while...\n")
  pca_analysis <- prcomp(listaintermedia,center=T,scale=T)
  if(is.null(annotation)){
    resultado <- list(pca_analysis,tratamientos,names(listadatos))
    names(resultado) <- c("pca","treatments","datasetnames")
    class(resultado) <- "pcaanalysis"
    return(resultado)
  } else {
    resultado <- list(pca_analysis,tratamientos,names(listadatos),annotation_table)
    names(resultado) <- c("pca","treatments","datasetnames","annotations")
    class(resultado) <- "pcaanalysis"
    return(resultado)
  }

}

pca_plot <- function(pca_list,treatment=1, compX=1, compY=2, compZ=NULL, plottype="Composite",useannot=FALSE,fontsizes=c(14,10,16,12)){
  library("plotly")
  ## Se ha puesto un tipo de plot por defecto porque si no, no sale nada, pero no devuelve error
  ## Se dara la opcion de escoger un tratamiento, que sera factor 1[1], factor 2 [2](si lo definio el usuario), interaccion[3], o combinar ambos[4].
  ## Falta corregir algunas funciones graficas, y temas de maquetado/leyendas, pero es cosa menor y no afecta a la funcionalidad
  ## Es guay, se puede llamar a una funcion de forma recursiva desde la misma funcion... (esto no lo sabia)

  ## Hay un bloque de vectores de "imagen" con los sets de colores y de simbolos. Como habra que cambiarlos, porque no seran
  ## todos del agrado de la mayoria, pues mejor que se definan fuera de todo.

  ## Plottype: Aqui estan, o estaran, todos los graficos que pueden ser utiles. Hay unos cuantos, la gente lega no los usara,
  ## pero es interesante mostrarlos.
  #             "Composite" <- Collage con los 4 plots mas representativos del analisis
  #             "Score" <- representa los individuos (1 punto por muestra), considerando todos los niveles. Unico que permitira 3D
  #             "Scree" <- Scree plot similar al PCA
  #             "Var" <- Muestrea las correlaciones de variables, permite escoger dataset y anotacion
  #             "Biplot" <- Biplots, a implementar
  #             "All" <- Devuelve una lista con todos los graficos anteriores.


  ## ********Esto hay que mejorarlo un poco controlando todos los niveles del objeto y el tipo de grafico*******
  ## Validamos que el objeto de verdad sea del tipo que tiene que ser:
  ##
  if(class(pca_list)!= "pcaanalysis")
    stop("PCA analysis object is expected, please run pca_analysis first")
  if(!inherits(pca_list[[1]], "prcomp"))
    stop("PCA analysis object seems to be wrong. Please re-run pca_analysis")
  if(plottype %in% c("Composite","Score","Scree","Var","Biplot","All")==FALSE)
    stop("Please select a valid plot type")
  if(useannot==T&&is.null(pca_list$annotations)){
    cat(paste("\nAnnotation matrix not provided. Considering useannot=F\n", "Please check/repeat pca_analysis() including an annotation matrix\n","Read help or follow the tutorial for more information.\n",sep=""))
    useannot=FALSE
  }
  #...

  ##
  ## Definimos los objetos necesarios para dibujar
  ##

  pca_object <- pca_list$pca
  treatmentstable <- pca_list$treatments
  annotations <- pca_list$annotations
  datasetnames <- pca_list$datasetnames
  treatmentnames <- colnames(treatmentstable)

  #Extraemos las varianzas explicadas. Esto lo necesitamos para las etiquetas de ejes, y para exportar.
  explainedvariance <- round((pca_object$sdev^2)*100/sum(pca_object$sdev^2),2)  #Sacamos la varianza explicada por cada componente

  #Plot sencillo, Scores
  if(plottype=="Score"){
    if(treatment!=4){  #Asi seleccionamos todas las visualizaciones que solo tengan un factor a mostrar.
      if(treatment==1) tratamientocolores <-as.vector(treatmentstable[,1])
      if(treatment==2) tratamientocolores <-as.vector(treatmentstable[,2])
      if(treatment==3) tratamientocolores <- as.vector(treatmentstable[,3])

      if(is.null(compZ)){ # Plot 2D, 1 factor.
        pca_score_plot <- plot_ly(as.data.frame(pca_object$x),
                                  x = ~pca_object$x[,compX],
                                  y = ~pca_object$x[,compY],
                                  hoverinfo ="text",
                                  text = paste("Sample:", rownames(pca_object$x), sep=" "),
                                  showlegend = TRUE,
                                  type="scatter",
                                  mode = "markers",
                                  color = factor(tratamientocolores,labels = c(unique(tratamientocolores))),
                                  colors = paleta_tratamientos[1:length(unique(tratamientocolores))],
                                  #symbol = factor(rep(c(1:nrow(mcoin_object$mcoa$lambda)),each=nrow(mcoin_object$mcoa$SynVar)), labels=datasetnames),
                                  symbols = simbolos[1:length(unique(tratamientocolores))],
                                  marker = list(size = 11)) %>%
          layout(xaxis = list(title = paste("Component",compX,"-",explainedvariance[compX],"% of variance", sep=" "),titlefont=list(size=fontsizes[1]),tickfont=list(size=fontsizes[2])),
                 yaxis = list(title = paste("Component",compY,"-",explainedvariance[compY],"% of variance",sep=" "),titlefont=list(size=fontsizes[1]),tickfont=list(size=fontsizes[2])),
                 title = "PCA ScorePlot",
                 titlefont = list(size=fontsizes[3]),
                 legend = list(font = list(size=fontsizes[4])))

        return(pca_score_plot)
      }

      if(!is.null(compZ)){ # Plot 3D, 1 factor.
        pca_score_plot <- plot_ly(as.data.frame(pca_object$x),
                                  x = ~pca_object$x[,compX],
                                  y = ~pca_object$x[,compY],
                                  z = ~pca_object$x[,compZ],
                                  hoverinfo ="text",
                                  text = paste("Sample:", rownames(pca_object$x), sep=" "),
                                  showlegend = TRUE,
                                  type="scatter3d",
                                  mode = "markers",
                                  color = factor(tratamientocolores,labels = c(unique(tratamientocolores))),
                                  colors = paleta_tratamientos[1:length(unique(tratamientocolores))],
                                  #symbol = factor(rep(c(1:nrow(mcoin_object$mcoa$lambda)),each=nrow(mcoin_object$mcoa$SynVar)), labels=datasetnames),
                                  symbols = simbolos[1:length(unique(tratamientocolores))],
                                  marker = list(size = 11)) %>%
          layout(scene=list(xaxis = list(title = paste("Component",compX,"-",explainedvariance[compX],"% of variance", sep=" "),titlefont=list(size=fontsizes[1]),tickfont=list(size=fontsizes[2])),
                 yaxis = list(title = paste("Component",compY,"-",explainedvariance[compY],"% of variance",sep=" "),titlefont=list(size=fontsizes[1]),tickfont=list(size=fontsizes[2])),
                 zaxis = list(title = paste("Component",compZ,"-",explainedvariance[compZ],"% of variance",sep=" "),titlefont=list(size=fontsizes[1]),tickfont=list(size=fontsizes[2]))),
                 title = "PCA ScorePlot",
                 titlefont = list(size=fontsizes[3]),
                 legend = list(font = list(size=fontsizes[4])))

        return(pca_score_plot)
      }

    }
    if(treatment==4){  #Asi seleccionamos las visualizaciones con 2 factores a mostrar

      if(is.null(compZ)){ # Plot 2D, 2 factores.
        pca_score_plot <- plot_ly(as.data.frame(pca_object$x),
                                  x = ~pca_object$x[,compX],
                                  y = ~pca_object$x[,compY],
                                  hoverinfo ="text",
                                  text = paste("Sample:", rownames(pca_object$x), sep=" "),
                                  showlegend = TRUE,
                                  type="scatter",
                                  mode = "markers",
                                  color = factor(treatmentstable[,1], labels = unique(treatmentstable[,1])),
                                  colors = paleta_tratamientos[1:length(unique(treatmentstable[,1]))],
                                  symbol = factor(treatmentstable[,2], labels=unique(treatmentstable[,2])),
                                  symbols = simbolos[1:length(unique(treatmentstable[,2]))],
                                  marker = list(size = 11)) %>%
          layout(xaxis = list(title = paste("Component",compX,"-",explainedvariance[compX],"% of variance", sep=" "),titlefont=list(size=fontsizes[1]),tickfont=list(size=fontsizes[2])),
                 yaxis = list(title = paste("Component",compY,"-",explainedvariance[compY],"% of variance",sep=" "),titlefont=list(size=fontsizes[1]),tickfont=list(size=fontsizes[2])),
                 title = "PCA ScorePlot",
                 titlefont = list(size=fontsizes[3]),
                 legend = list(font = list(size=fontsizes[4])))

        return(pca_score_plot)


      }

      if(!is.null(compZ)){ # Plot 3D, 2 factores.
        pca_score_plot <- plot_ly(as.data.frame(pca_object$x),
                                  x = ~pca_object$x[,compX],
                                  y = ~pca_object$x[,compY],
                                  z = ~pca_object$x[,compZ],
                                  hoverinfo ="text",
                                  text = paste("Sample:", rownames(pca_object$x), sep=" "),
                                  showlegend = TRUE,
                                  type="scatter3d",
                                  mode = "markers",
                                  color = factor(treatmentstable[,1], labels = unique(treatmentstable[,1])),
                                  colors = paleta_tratamientos[1:length(unique(treatmentstable[,1]))],
                                  symbol = factor(treatmentstable[,2], labels=unique(treatmentstable[,2])),
                                  symbols = simbolos[1:length(unique(treatmentstable[,2]))],
                                  marker = list(size = 11)) %>%
          layout(scene=list(xaxis = list(title = paste("Component",compX,"-",explainedvariance[compX],"% of variance", sep=" "),titlefont=list(size=fontsizes[1]),tickfont=list(size=fontsizes[2])),
                 yaxis = list(title = paste("Component",compY,"-",explainedvariance[compY],"% of variance",sep=" "),titlefont=list(size=fontsizes[1]),tickfont=list(size=fontsizes[2])),
                 zaxis = list(title = paste("Component",compZ,"-",explainedvariance[compZ],"% of variance",sep=" "),titlefont=list(size=fontsizes[1]),tickfont=list(size=fontsizes[2]))),
                 title = "PCA ScorePlot",
                 titlefont = list(size=fontsizes[3]),
                 legend = list(font = list(size=fontsizes[4])))

        return(pca_score_plot)

      }
    }
  }

  if(plottype=="Scree"){
    explainedvariance <- pca_object$sdev^2  #Sacamos la varianza explicada por cada tratamiento
    screeplotdata <- data.frame(PCs=paste("PC",c(1:length(explainedvariance)),sep=""), ExplVar=round(explainedvariance*100/sum(explainedvariance),2), Cumulat=round(cumsum(explainedvariance)*100/sum(explainedvariance),2),stringsAsFactors = FALSE) #Creamos tabla de datos con % varianza y % var acumulado
    screeplotdata$PCs <- factor(screeplotdata$PCs, levels = unique(screeplotdata$PCs)[order(screeplotdata$ExplVar, decreasing = TRUE)]) # Ordenamos los PCs segun valor

    pca_screeplot <- plot_ly(screeplotdata) %>%
      add_trace(x=~PCs,y=~ExplVar,type="bar",name="% Explained Variance",
                marker = list(color = '#C9EFF9')) %>%
      add_trace(x=~PCs, y=~Cumulat,type="scatter",mode="lines+markers",yaxis = "y2",name="Cumulative Proportion",
                line = list(color = '#45171D')) %>%
      layout(title = "PCA Screeplot",
             titlefont = list(size=fontsizes[3]),
             xaxis = list(title = "",titlefont=list(size=fontsizes[1]),tickfont=list(size=fontsizes[2])),
             yaxis = list(side = 'left', title = "% of Explained Variance", showgrid = F, zeroline = TRUE,titlefont=list(size=fontsizes[1]),tickfont=list(size=fontsizes[2])),
             yaxis2 = list(side = 'right', overlaying = "y", title = "Cumulative", showgrid = FALSE, zeroline = FALSE,titlefont=list(size=fontsizes[1]),tickfont=list(size=fontsizes[2])),
             legend = list(x = 0.7, y = 0.1),font = list(size=fontsizes[4]))

    return(pca_screeplot)

  }

  if(plottype=="Var"){
    #Generamos anotaciones
    annotations2 <- as.data.frame(row.names(pca_object$rotation))
    colnames(annotations2)<-"IDENTIFIER"
    annotations <- left_join(annotations2,as.data.frame(annotations),by="IDENTIFIER") #con esto asignamos cada anotación a la variable correspondiente por rowname. Asi no hay fallo.
    annotations <- as.matrix(annotations) #Por alguna razon, sin esto, luego no se ordenan los factores ni los colores.
    rm(annotations2)

    if(useannot==FALSE){
      pca_VarPlot <- plot_ly(as.data.frame(pca_object$rotation),
                             x = ~pca_object$rotation[,compX],
                             y = ~pca_object$rotation[,compY],
                             text = paste("ID:", rownames(pca_object$rotation)),
                             type="scatter",
                             mode = "markers",
                             marker = list(size = 6, opacity=1),
                             colors =paleta_niveles[1:length(datasetnames)])

        pca_VarPlot <- layout(pca_VarPlot,
                              xaxis = list(title = paste("Component",compX,"-",explainedvariance[compX],"% of variance", sep=" "),titlefont=list(size=fontsizes[1]),tickfont=list(size=fontsizes[2])),
                              yaxis = list(title = paste("Component",compY,"-",explainedvariance[compY],"% of variance",sep=" "),titlefont=list(size=fontsizes[1]),tickfont=list(size=fontsizes[2])),
                              title = "PCA Variable Plot",
                              titlefont = list(size=fontsizes[3]),
                              legend = list(font = list(size=fontsizes[4])))


      return(pca_VarPlot)
    }
    if(useannot==TRUE){
      #chapuza para las anotaciones. eliminamos duplicados pa porsi. en base al ID
      annotations<-annotations[!duplicated(annotations[,1]),]
      
      #Sacamos los niveles del mapman, y las etiquetas
      mapmanlevels <- as.numeric(unique(annotations[order(as.numeric(annotations[,3])),3]))
      mapmannames <- annotations[order(as.numeric(annotations[,3])),]
      mapmannames <- as.vector(unique(mapmannames[,4]))

      pca_VarPlot <- plot_ly(as.data.frame(pca_object$rotation),
                             x = ~pca_object$rotation[,compX],
                             y = ~pca_object$rotation[,compY],
                             text = paste("ID:", rownames(pca_object$rotation),"</br> Description: ",annotations[,2],"</br> Mercator Bin: ",annotations[,4]),
                             type="scatter",
                             mode = "markers",
                             marker = list(size = 6, opacity=1),
                             color = factor(annotations[,4], levels=mapmannames),
                             colors = paleta_mapman[mapmanlevels])

      pca_VarPlot <- layout(pca_VarPlot,
                           xaxis = list(title = paste("Component",compX,"-",explainedvariance[compX],"% of variance", sep=" "),titlefont=list(size=fontsizes[1]),tickfont=list(size=fontsizes[2])),
                           yaxis = list(title = paste("Component",compY,"-",explainedvariance[compY],"% of variance",sep=" "),titlefont=list(size=fontsizes[1]),tickfont=list(size=fontsizes[2])),
                           title = "PCA Variable Plot",
                           titlefont = list(size=fontsizes[3]),
                           legend = list(font = list(size=fontsizes[4])))

      return(pca_VarPlot)
    }
  }

  if(plottype=="Biplot"){
    #Generamos anotaciones
    annotations2 <- as.data.frame(row.names(pca_object$rotation))
    colnames(annotations2)<-"IDENTIFIER"
    annotations <- left_join(annotations2,as.data.frame(annotations),by="IDENTIFIER") #con esto asignamos cada anotación a la variable correspondiente por rowname. Asi no hay fallo.
    annotations <- as.matrix(annotations) #Por alguna razon, sin esto, luego no se ordenan los factores ni los colores.
    rm(annotations2)

    ##### Rollo biplot: Aquí deberemos escalar y ajustar varianzas tanto de loadigns como de scores. Hay forma canónica de hacerlo, y es la siguiente
    scores <- pca_object$x  # sacamos la matriz con los scores
    scale <- 0.9 #factor escala de los loadings, esto es importante. En principio se deja así, pero igual hay que probar visualizaciones en rangos 0.8-1. El biplot básico de R toma 0.8
    lam <- pca_object$sdev[c(compX,compY)] # Raiz cuadrada de evalores de los componentes seleccionados
    n <- nrow(pca_object$x) # numero de filas de los escores
    lam <- lam *sqrt(n) # Escalado de lam, leer debajo
    # Los scores de prcomp tienen varianza igual al autovalor. Si los dividimos por la raiz cuadrada del autovalor escalaremos por unidad de varianza.
    # Si queremos escalar por la suma de cuadrados, hay que multiplicarlo por la raiz de n

    x <- scores[,c(compX,compY)]/lam  # Escalamos los scores de las muestras.
    y <- pca_object$rotation[,c(compX,compY)]*lam # Escalamos las cargas.

    unsigned.range <- function(x) c(-abs(min(x, na.rm = TRUE)), abs(max(x, na.rm = TRUE))) # sacamos cuales son los rangos para cada eje de scores y loadings
    rangx1 <- unsigned.range(x[, 1L])
    rangx2 <- unsigned.range(x[, 2L])
    rangy1 <- unsigned.range(y[, 1L])
    rangy2 <- unsigned.range(y[, 2L])

    ratioejex <- rangx1/rangy1
    ratioejey <- rangx2/rangy2

    ratio <- max(ratioejex,ratioejey)
    y=y*ratio*scale
    ### Fin de las normalizaciones para el biplot


    if(treatment!=4){  #Asi seleccionamos todas las visualizaciones que solo tengan un factor a mostrar.
      if(treatment==1) tratamientocolores <-treatmentstable[,1]
      if(treatment==2) tratamientocolores <-treatmentstable[,2]
      if(treatment==3) tratamientocolores <-treatmentstable[,3]

      if(useannot==FALSE){
        #Sacamos los niveles del mapman, y las etiquetas
        paletaintegrada <- c("#466fc7", paleta_tratamientos[1:length(unique(tratamientocolores))])

        pca_biplot <- plot_ly()
        pca_biplot <- add_markers(pca_biplot,
                                  x = ~y[,1],
                                  y = ~y[,2],
                                  text = paste("ID:", rownames(pca_object$rotation)),
                                  type="scatter",
                                  mode = "markers",
                                  marker = list(size = 6, opacity=1),
                                  color = as.factor("Loadings"),
                                  colors = paletaintegrada)
        pca_biplot <- add_markers(pca_biplot,
                                  x = ~x[,1],
                                  y = ~x[,2],
                                  inherit = FALSE,
                                  hoverinfo ="text",
                                  text = paste("Sample:", rownames(x), sep=" "),
                                  showlegend = TRUE,
                                  type="scatter",
                                  mode = "markers",
                                  color = factor(tratamientocolores,labels = c(unique(tratamientocolores))),
                                  colors = paletaintegrada[2:length(paletaintegrada)],
                                  marker = list(size = 11),
                                  yaxis="y",
                                  xaxis="x")

          pca_biplot <- layout(pca_biplot,
                               xaxis = list(title = paste("Component",compX,"-",explainedvariance[compX],"% of variance", sep=" "),titlefont=list(size=fontsizes[1]),tickfont=list(size=fontsizes[2])),
                               yaxis = list(title = paste("Component",compY,"-",explainedvariance[compY],"% of variance",sep=" "),titlefont=list(size=fontsizes[1]),tickfont=list(size=fontsizes[2])),
                               title = "PCA Biplot",
                               titlefont = list(size=fontsizes[3]),
                               legend = list(font = list(size=fontsizes[4])))

        return(pca_biplot)
      }
      if(useannot==TRUE){
        #chapuza para las anotaciones. eliminamos duplicados pa porsi. en base al ID
        annotations<-annotations[!duplicated(annotations[,1]),]
        
        #Sacamos los niveles del mapman, y las etiquetas
        mapmanlevels <- as.numeric(unique(annotations[order(as.numeric(annotations[,3])),3]))
        mapmannames <- annotations[order(as.numeric(annotations[,3])),]
        mapmannames <- as.vector(unique(mapmannames[,4]))
        paletaintegrada <- c(paleta_mapman[mapmanlevels], paleta_tratamientos[1:length(unique(tratamientocolores))])

        pca_biplot <- plot_ly()
          pca_biplot <- add_markers(pca_biplot,
                      x = ~y[,1],
                      y = ~y[,2],
                      text = paste("ID:", rownames(pca_object$rotation),"</br> Description: ",annotations[,2],"</br> Mercator Bin: ",annotations[,4]),
                      type="scatter",
                      mode = "markers",
                      marker = list(size = 6, opacity=1),
                      color = factor(annotations[,4], levels=mapmannames),
                      colors = paletaintegrada)
            pca_biplot <- add_markers(pca_biplot,
                      x = ~x[,1],
                      y = ~x[,2],
                      inherit = FALSE,
                      hoverinfo ="text",
                      text = paste("Sample:", rownames(x), sep=" "),
                      showlegend = TRUE,
                      type="scatter",
                      mode = "markers",
                      color = factor(tratamientocolores,labels = c(unique(tratamientocolores))),
                      colors = paletaintegrada[(length(paletaintegrada)-length(unique(tratamientocolores))):(length(paletaintegrada))],
                      marker = list(size = 11),
                      yaxis="y",
                      xaxis="x")

          pca_biplot <- layout(pca_biplot,
                 xaxis = list(title = paste("Component",compX,"-",explainedvariance[compX],"% of variance", sep=" "),titlefont=list(size=fontsizes[1]),tickfont=list(size=fontsizes[2])),
                 yaxis = list(title = paste("Component",compY,"-",explainedvariance[compY],"% of variance",sep=" "),titlefont=list(size=fontsizes[1]),tickfont=list(size=fontsizes[2])),
                 title = "PCA Biplot",
                 titlefont = list(size=fontsizes[3]),
                 legend = list(font = list(size=fontsizes[4])))

        return(pca_biplot)
      }
    }

    if(treatment==4){  #Asi seleccionamos las visualizaciones con 2 factores a mostrar
      if(useannot==FALSE){
        #Sacamos los niveles del mapman, y las etiquetas
        paletaintegrada <- c("#466fc7", paleta_tratamientos[1:length(unique(treatmentstable[,1]))])

        pca_biplot <- plot_ly()
        pca_biplot <- add_markers(pca_biplot,
                                  x = ~y[,1],
                                  y = ~y[,2],
                                  text = paste("ID:", rownames(pca_object$rotation)),
                                  type="scatter",
                                  mode = "markers",
                                  marker = list(size = 6, opacity=1),
                                  color = as.factor("Loadings"),
                                  colors = paletaintegrada)

        pca_biplot <- add_trace(pca_biplot,
                              x = ~x[,1],
                  y = ~x[,2],
                  hoverinfo ="text",
                  text = paste("Sample:", rownames(x), sep=" "),
                  showlegend = TRUE,
                  type="scatter",
                  mode = "markers",
                  color = factor(treatmentstable[,1], labels = unique(treatmentstable[,1])),
                  colors = paletaintegrada[2:length(paletaintegrada)],
                  symbol = factor(treatmentstable[,2], labels=unique(treatmentstable[,2])),
                  symbols = simbolos[1:length(unique(treatmentstable[,2]))],
                  marker = list(size = 11),
                  inherit = F)

        pca_biplot <- layout(pca_biplot,
                             xaxis = list(title = paste("Component",compX,"-",explainedvariance[compX],"% of variance", sep=" "),titlefont=list(size=fontsizes[1]),tickfont=list(size=fontsizes[2])),
                             yaxis = list(title = paste("Component",compY,"-",explainedvariance[compY],"% of variance",sep=" "),titlefont=list(size=fontsizes[1]),tickfont=list(size=fontsizes[2])),
                             title = "PCA Biplot",
                             titlefont = list(size=fontsizes[3]),
                             legend = list(font = list(size=fontsizes[4])))

        return(pca_biplot)

      }
      if(useannot==TRUE){
        #chapuza para las anotaciones. eliminamos duplicados pa porsi. en base al ID
        annotations<-annotations[!duplicated(annotations[,1]),]
        
        #Sacamos los niveles del mapman, y las etiquetas
        mapmanlevels <- as.numeric(unique(annotations[order(as.numeric(annotations[,3])),3]))
        mapmannames <- annotations[order(as.numeric(annotations[,3])),]
        mapmannames <- as.vector(unique(mapmannames[,4]))
        paletaintegrada <- c(paleta_mapman[mapmanlevels], paleta_tratamientos[1:length(unique(treatmentstable[,1]))])

        pca_biplot <- plot_ly()
        pca_biplot <- add_markers(pca_biplot,
                                  x = ~y[,1],
                                  y = ~y[,2],
                                  text = paste("ID:", rownames(pca_object$rotation),"</br> Description: ",annotations[,2],"</br> Mercator Bin: ",annotations[,4]),
                                  type="scatter",
                                  mode = "markers",
                                  marker = list(size = 6, opacity=1),
                                  color = factor(annotations[,4], levels=mapmannames),
                                  colors = paletaintegrada)
        pca_biplot <- add_markers(pca_biplot,
                                  x = ~x[,1],
                                  y = ~x[,2],
                                  inherit = FALSE,
                                  hoverinfo ="text",
                                  text = paste("Sample:", rownames(x), sep=" "),
                                  showlegend = TRUE,
                                  type="scatter",
                                  mode = "markers",
                                  color = factor(treatmentstable[,1],labels = c(unique(treatmentstable[,1]))),
                                  colors = paletaintegrada[(length(paletaintegrada)-length(unique(treatmentstable[,1]))):(length(paletaintegrada))],
                                  marker = list(size = 11),
                                  yaxis="y",
                                  xaxis="x",
                                  symbol = factor(treatmentstable[,2], labels=unique(treatmentstable[,2])),
                                  symbols = simbolos[1:length(unique(treatmentstable[,2]))],
                                  marker = list(size = 11))

        pca_biplot <- layout(pca_biplot,
                             xaxis = list(title = paste("Component",compX,"-",explainedvariance[compX],"% of variance", sep=" "),titlefont=list(size=fontsizes[1]),tickfont=list(size=fontsizes[2])),
                             yaxis = list(title = paste("Component",compY,"-",explainedvariance[compY],"% of variance",sep=" "),titlefont=list(size=fontsizes[1]),tickfont=list(size=fontsizes[2])),
                             title = "PCA Biplot",
                             titlefont = list(size=fontsizes[3]),
                             legend = list(font = list(size=fontsizes[4])))


        return(pca_biplot)
      }
    }

  }

  if(plottype=="Composite"){
    if(!is.null(compZ)) stop("Composite plot can only draw 2D plots. compZ must be null")
    #Representamos los eigenvalues, y luego como sube la proporcion acumulada.
    p1 <- pca_plot(pca_list,treatment, compX, compY, plottype = "Score", fontsizes = fontsizes)
    p2 <- pca_plot(pca_list,treatment, compX, compY, plottype = "Var", useannot=FALSE,fontsizes = fontsizes)
    p3 <- pca_plot(pca_list,treatment, compX, compY, plottype = "Biplot",useannot = useannot,fontsizes = fontsizes)
    #p4 <- pca_plot(pca_list,treatment, compX, compY, plottype = "Scree")
    #mcoa_ScorePlot <- subplot(subplot(p1,p2),subplot(p3,p4),nrows=2,shareX = F,shareY = F,margin = 0.04)
    pca_CompositeplotPlot <- subplot(subplot(p1,p2),subplot(p3),nrows=2,shareX = F,shareY = F,margin = 0.04)
    pca_CompositeplotPlot <- layout(pca_CompositeplotPlot, title = "PCA Analysis",titlefont = list(size=fontsizes[3]))
    return(pca_CompositeplotPlot)
  }

  if(plottype=="All"){
    #Aqui habra que crear una llamada recursiva a todos los graficos posibles y devolver una lista.
    #tambien recursivo pa los distintos niveles.
    if(!is.null(compZ)) stop("Please select only two dimensions for plotting all figures")
    cat("EXPORTING PCA ANALYSIS PLOTS\nA single core will be used. It may take a while...\n ")
    cat("Generating and saving Score plot. Step 1/5\n")
    export_plot(pca_plot(pca_list,treatment=treatment, compX=compX, compY=compY, plottype="Score",useannot=useannot),"pca_score.pdf")
    cat("Generating and saving Scree plot. Step 2/5\n")
    export_plot(pca_plot(pca_list,treatment=treatment, compX=compX, compY=compY, plottype="Scree",useannot=useannot),"pca_scree.pdf")
    cat("Generating and saving Var plot. Step 3/5\n")
    export_plot(pca_plot(pca_list,treatment=treatment, compX=compX, compY=compY, plottype="Var",useannot=useannot),"pca_var.pdf")
    cat("Generating and saving Biplot. Step 4/5\n")
    export_plot(pca_plot(pca_list,treatment=treatment, compX=compX, compY=compY, plottype="Biplot",useannot=useannot),"pca_biplot.pdf")
    cat("Generating and saving Composite plot. Step 5/5\n")
    export_plot(pca_plot(pca_list,treatment=treatment, compX=compX, compY=compY, plottype="Composite",useannot=useannot),"pca_composite.pdf")
  }
}
