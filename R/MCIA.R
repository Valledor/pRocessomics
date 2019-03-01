# ANALISIS DE COINERCIA MULTIPLE
#
# MCOA, require omicade4, ade4, plotly

# Requiere una lista con tantos niveles ?micos como vayamos a estudiar.
# Los individuos se ordenan en columnas, las variables en filas.
# Los nombres de todas las columnas han de ser iguales en los distintos elementos de la lista.
# y poco m?s. Ver ejemplo:

# Creamos lista con los datos ?micos antes de empezar. Se supone que ya la tendremos.
# Los elementos de la lista debern nombrarse de acuerdo a como el usuario haya nombrando los datasets

# nota hay que añadir algo para controlar los nas cuando hay puntos en los nombres de variables en la anotacion

##
## ANALISIS ESTADÍSTICO
##

mcia_analysis <- function(listadatos, annotation=NULL, initialrow=1, initialcolumn=3,treatment1col=1, treatment2col=2, treatment=1, omiclevel="all", keptaxes=5){
  require("ade4")
  require("omicade4")
  require("plotly")
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
  if(treatment1col==treatment2col) tratamientos <-cbind(listadatos[[1]][initialrow:nrow(listadatos[[1]]),treatment1col],listadatos[[1]][initialrow:nrow(listadatos[[1]]),treatment1col])
  if(treatment1col!=treatment2col) tratamientos <-cbind(listadatos[[1]][initialrow:nrow(listadatos[[1]]),treatment1col],listadatos[[1]][initialrow:nrow(listadatos[[1]]),treatment2col])
  tratamientos <- cbind(tratamientos, paste(listadatos[[1]][initialrow:nrow(listadatos[[1]]),treatment1col],listadatos[[1]][initialrow:nrow(listadatos[[1]]),treatment2col], sep="&"))
  
  if(is.null(omiclevel)==TRUE) omiclevel="all"
  listadatos <- lapply(listadatos,function(x) x[initialrow:nrow(x),initialcolumn:ncol(x)])
  
  # Seleccionamos los niveles omicos datos
  if(omiclevel=="all") {
    listadatos2<- listadatos
    listadatos2<-lapply(listadatos2, function(x) RemoveEmptyColumns(x,1,1))#Esto es para evitar cosas raras.
  }  
  if(omiclevel!="all") {
    sapply(omiclevel, function(x) if(x %in% datasetnames==FALSE) stop("Error: Please select a valid datasetname"))
    listadatos2<-listadatos[omiclevel]
  }
  
  #Extraemos los valores significativos de las tablas
  listadatos2 <- lapply(listadatos2, function(x) x[initialrow:nrow(x),initialcolumn:ncol(x)])
  variablenames <- lapply(listadatos2, function(x) colnames(x))
  variablenames <-unlist(variablenames)
  

  if(!is.null(annotation)){ #Creamos una tabla de anotaciones, a nuestra manera, para plotear luego. Chapuza, pero funciona.
    annotation_table <- annotation_mapman(listadatos2,annotation)
  }
  
  # Ojo que para este analisis las variables van en filas.
  listadatos2<-lapply(listadatos2, function(x) t(x))
  
  # Comprobamos que haya al menos dos niveles, todos los elementos tengan las mismas columnas, y que todas se llamen igual.
  if(typeof(listadatos2 ) != "list")
    stop("A list of matrices or dataframes corresponding to each level is expected")
  if(is.null(names(listadatos2)))
    stop("Names of elements in list are missing")
  if(length(listadatos2 ) <= 1)
    stop("At least two levels are necessary to perform this test")
  #if(all(sapply(listadatos2, dim)[2,]==sapply(listadatos2, dim)[2,1])==FALSE)
    #stop("The different matrices have an unequal number of individuals")
  #if(all(apply((x <- sapply(listadatos2, colnames))[,-1], 2, function(y) identical(y, x[,1])))==FALSE)
    #stop("Individuals have different names across levels. Please check original matrices")
  
  
  # Hacemos el análisis 
  mcoin_analysis <- mcia(listadatos2, cia.nf=keptaxes, cia.scan=F) 
  
  # Retornamos una lista con dos objetos, en análisis, y las 2 columntas de la tabla con los tratamientos para plotear luego
  if(is.null(annotation)){
    resultado <- list(mcoin_analysis,tratamientos,names(listadatos))
    names(resultado) <- c("mcoa","treatments","datasetnames")
    class(resultado) <- "mcoaanalysis"
  return(resultado)
  } else {
    resultado <- list(mcoin_analysis,tratamientos,names(listadatos),annotation_table)
    names(resultado) <- c("mcoa","treatments","datasetnames","annotations")
    class(resultado) <- "mcoaanalysis"
    return(resultado)
  }
}


##
## Funcion para generar todos los tipos de dibujos posibles.
##

mcia_plot <- function(mcoin_list,treatment=1, compX=1, compY=2, compZ=NULL, plottype="Composite",useannot=FALSE,levelchoice="all",fontsizes=c(14,10,16,12)){ 
  ## Se ha puesto un tipo de plot por defecto porque si no, no sale nada, pero no devuelve error
  ## Se dara la opcion de escoger un tratamiento, que sera factor 1[1], factor 2 [2](si lo definio el usuario), interaccion[3], o combinar ambos[4].
  ## Falta corregir algunas funciones graficas, y temas de maquetado/leyendas, pero es cosa menor y no afecta a la funcionalidad
  ## Es guay, se puede llamar a una funcion de forma recursiva desde la misma funcion... (esto no lo sabia)
  
  ## Hay un bloque de vectores de "imagen" con los sets de colores y de simbolos. Como habra que cambiarlos, porque no seran
  ## todos del agrado de la mayoria, pues mejor que se definan fuera de todo.
  
  ## Plottype: Aqui estan, o estaran, todos los graficos que pueden ser utiles. Hay unos cuantos, la gente lega no los usara, 
  ## pero es interesante mostrarlos.
  #             "Composite" <- Collage con los 4 plots mas representativos del analisis
  #             "SynVar" <- representa los individuos (1 punto por muestra), considerando todos los niveles. Unico que permitira 3D
  #             "Scree" <- Scree plot similar al PCA
  #             "Var" <- Muestrea las correlaciones de variables, permite escoger dataset y anotacion
  #             "CoIn" <- Representacion del analisis de coinercia propiamente dicho
  #             "Biplot" <- Listo,pero con fallos de paletas
  #             "Pseudoeigen" <- Muestra lo correlacionadas o no que estan los datasets
  #             "All" <- Devuelve una lista con todos los graficos anteriores.
  
  
  ## ********Esto hay que mejorarlo un poco controlando todos los niveles del objeto y el tipo de grafico*******
  ## Validamos que el objeto de verdad sea del tipo que tiene que ser:
  ##
  if(class(mcoin_list)!= "mcoaanalysis")
    stop("MCOA analysis object is expected, please run pca_analysis first")
  if(!inherits(mcoin_list[[1]], "mcia"))
    stop("mcia object expected, please run mcoin_analysis first")
  if(plottype %in% c("Composite","SynVar","Scree","Var", "VarAnnot","CoIn","Pseudoeigen","All","Biplot")==FALSE)
    stop("Please select a valid plot type")
  if(useannot==T&&is.null(mcoin_list$annotations)){
    cat(paste("\nAnnotation matrix not provided. Considering useannot=F\n", "Please check/repeat mcoin_analysis() including an annotation matrix\n","Read help or follow the tutorial for more information.\n",sep=""))
    useannot=FALSE
  }
  #...

  
  ##
  ## Definimos los objetos necesarios para dibujar
  ##
  
  mcoin_object <- mcoin_list$mcoa
  treatmentstable <- mcoin_list$treatments
  annotations <- mcoin_list$annotations
  datasetnames <- mcoin_list$datasetnames
  treatmentnames <- colnames(treatmentstable)
  
  # Determinamos el n?mero de niveles que tenemos, para no depender de la lisa empleada en mcoin_analysis, lo sacamos de los propios resultados
  nivelesomicos <- nrow(mcoin_object$mcoa$lambda) #En teoria debemos pasarle los tratamientos...
  variablesxnivel <- as.vector(mcoin_object$mcoa$blo) #Aquí sabemos cuantas variables tenemos por nivel
  niveles <- length(variablesxnivel)
  
  #Plot sencillo, SynVar
  if(plottype=="SynVar"){
                          
    if(treatment!=4){  #Asi seleccionamos todas las visualizaciones que solo tengan un factor a mostrar.
      if(treatment==1) tratamientocolores <-as.vector(treatmentstable[,1])
      if(treatment==2) tratamientocolores <-as.vector(treatmentstable[,2])
      if(treatment==3) tratamientocolores <- as.vector(treatmentstable[,3])
      
      if(is.null(compZ)){ # Plot 2D, 1 factor.
        mcoa_ScorePlot <- plot_ly(as.data.frame(mcoin_object$mcoa$SynVar), 
                                  x = ~mcoin_object$mcoa$SynVar[,compX], 
                                  y = ~mcoin_object$mcoa$SynVar[,compY],
                                  hoverinfo ="text",
                                  text = paste("Sample:", rownames(mcoin_object$mcoa$SynVar), sep=" "),
                                  showlegend = TRUE,
                                  type="scatter",
                                  mode = "markers", 
                                  color = factor(tratamientocolores,labels = c(unique(tratamientocolores))),
                                  colors = paleta_tratamientos[1:length(unique(tratamientocolores))],
                                  #symbol = factor(rep(c(1:nrow(mcoin_object$mcoa$lambda)),each=nrow(mcoin_object$mcoa$SynVar)), labels=datasetnames),
                                  symbols = simbolos[1:length(unique(tratamientocolores))], 
                                  marker = list(size = 11)) %>%
          layout(xaxis = list(title = paste("Component",compX,sep=" "),titlefont=list(size=fontsizes[1]),tickfont=list(size=fontsizes[2])),
                 yaxis = list(title = paste("Component",compY,sep=" "),titlefont=list(size=fontsizes[1]),tickfont=list(size=fontsizes[2])),
                 title = "MCIA Synthetic Scores Plot",
                 titlefont = list(size=fontsizes[3]),
                 legend = list(font = list(size=fontsizes[4])))
  
        return(mcoa_ScorePlot)
        }
     
      if(!is.null(compZ)){ # Plot 3D, 1 factor.
        mcoa_ScorePlot <- plot_ly(as.data.frame(mcoin_object$mcoa$SynVar), 
                                  x = ~mcoin_object$mcoa$SynVar[,compX], 
                                  y = ~mcoin_object$mcoa$SynVar[,compY],
                                  z = ~mcoin_object$mcoa$SynVar[,compZ],
                                  hoverinfo ="text",
                                  text = paste("Sample:", rownames(mcoin_object$mcoa$SynVar), sep=" "),
                                  showlegend = TRUE,
                                  type="scatter3d",
                                  mode = "markers", 
                                  color = factor(tratamientocolores,labels = c(unique(tratamientocolores))),
                                  colors = paleta_tratamientos[1:length(unique(tratamientocolores))],
                                  #symbol = factor(rep(c(1:nrow(mcoin_object$mcoa$lambda)),each=nrow(mcoin_object$mcoa$SynVar)), labels=datasetnames),
                                  symbols = simbolos[1:length(unique(tratamientocolores))], 
                                  marker = list(size = 11)) %>%

          layout(scene=list(xaxis = list(title = paste("Component",compX,sep=" "),titlefont=list(size=fontsizes[1]),tickfont=list(size=fontsizes[2])),
                            yaxis = list(title = paste("Component",compY,sep=" "),titlefont=list(size=fontsizes[1]),tickfont=list(size=fontsizes[2])),
                            zaxis = list(title = paste("Component",compZ,sep=" "),titlefont=list(size=fontsizes[1]),tickfont=list(size=fontsizes[2]))),
                 title = "MCIA Synthetic Scores Plot",
                 titlefont = list(size=fontsizes[3]),
                 legend = list(font = list(size=fontsizes[4])))
        
        return(mcoa_ScorePlot)
      }
       
    }
    if(treatment==4){  #Asi seleccionamos las visualizaciones con 2 factores a mostrar
      
      if(is.null(compZ)){ # Plot 2D, 2 factores.
        mcoa_ScorePlot <- plot_ly(as.data.frame(mcoin_object$mcoa$SynVar), 
                                  x = ~mcoin_object$mcoa$SynVar[,compX], 
                                  y = ~mcoin_object$mcoa$SynVar[,compY],
                                  hoverinfo ="text",
                                  text = paste("Sample:", rownames(mcoin_object$mcoa$SynVar), sep=" "),
                                  showlegend = TRUE,
                                  type="scatter",
                                  mode = "markers", 
                                  color = factor(treatmentstable[,1], labels = unique(treatmentstable[,1])),
                                  colors = paleta_tratamientos[1:length(unique(treatmentstable[,1]))],
                                  symbol = factor(treatmentstable[,2], labels=unique(treatmentstable[,2])),
                                  symbols = simbolos[1:length(unique(treatmentstable[,2]))], 
                                  marker = list(size = 11)) %>%
          layout(xaxis = list(title = paste("Component",compX,sep=" "),titlefont=list(size=fontsizes[1]),tickfont=list(size=fontsizes[2])),
                 yaxis = list(title = paste("Component",compY,sep=" "),titlefont=list(size=fontsizes[1]),tickfont=list(size=fontsizes[2])),
                 title = "MCIA Synthetic Scores Plot",
                 titlefont = list(size=fontsizes[3]),
                 legend = list(font = list(size=fontsizes[4])))
        
        return(mcoa_ScorePlot)
      }
      
      if(!is.null(compZ)){ # Plot 3D, 2 factores.
        mcoa_ScorePlot <- plot_ly(as.data.frame(mcoin_object$mcoa$SynVar), 
                                  x = ~mcoin_object$mcoa$SynVar[,compX], 
                                  y = ~mcoin_object$mcoa$SynVar[,compY],
                                  z = ~mcoin_object$mcoa$SynVar[,compZ],
                                  hoverinfo ="text",
                                  text = paste("Sample:", rownames(mcoin_object$mcoa$SynVar), sep=" "),
                                  showlegend = TRUE,
                                  type="scatter3d",
                                  mode = "markers", 
                                  color = factor(treatmentstable[,1], labels = unique(treatmentstable[,1])),
                                  colors = paleta_tratamientos[1:length(unique(treatmentstable[,1]))],
                                  symbol = factor(treatmentstable[,2], labels=unique(treatmentstable[,2])),
                                  symbols = simbolos[1:length(unique(treatmentstable[,2]))], 
                                  marker = list(size = 11)) %>%
          layout(scene=list(xaxis = list(title = paste("Component",compX,sep=" "),titlefont=list(size=fontsizes[1]),tickfont=list(size=fontsizes[2])),
                            yaxis = list(title = paste("Component",compY,sep=" "),titlefont=list(size=fontsizes[1]),tickfont=list(size=fontsizes[2])),
                            zaxis = list(title = paste("Component",compZ,sep=" "),titlefont=list(size=fontsizes[1]),tickfont=list(size=fontsizes[2]))),
                            title = "MCIA Synthetic Scores Plot",
                 titlefont = list(size=fontsizes[3]),
                 legend = list(font = list(size=fontsizes[4])))
        
        return(mcoa_ScorePlot)
      }
    }
  }
  
  if(plottype=="CoIn"){
    ## A partir de los scores sinteticos (la combinacion de los scores para los distintos niveles) vamos a definir unos segmentos
    ## (x1,y1),(x2,y2) para dibujar las conexiones de los distintos niveles omicos de cada muestra. Hay un valor fijo, el central,
    ## que es este SynVar, y luego tendremos tantos extremos como niveles omicos haya en su respectiva coordenada.
    
    # Definimos el tratamiento a usar para colorear
    if(treatment==1) tratamientocolores <-rep(treatmentstable[,1],nivelesomicos)
    if(treatment==2) tratamientocolores <-rep(treatmentstable[,2],nivelesomicos)
    if(treatment==3) tratamientocolores <- rep(paste(treatmentstable[,1],treatmentstable[,2],sep="&"),nivelesomicos)
    if(treatment==4) stop("CoIn plot do not accept two factors")
      #{tratamientocolores <- cbind(treatmentstable,paste(treatmentstable[,1],treatmentstable[,2],sep="-"))
      #colnames(tratamientocolores)[3] <- paste(treatmentnames[1],treatmentnames[2],sep="&")}
    
    segmentos <- cbind(mcoin_object$mcoa$SynVar[,compX],mcoin_object$mcoa$SynVar[,compY]) # Punto central X,Y para la combinacion de los componentes elegidos:
    segmentos_rep <- segmentos
    # Ahora la idea es crear una tabla de coordenadas en las que el valor centralha de repetirse tantas veces como niveles omicos tengamos:
    for(i in 2:nivelesomicos) segmentos_rep <- rbind(segmentos_rep, segmentos)
    # Añadimos las columnas con los scores de los distintos niveles ómicos. El orden se mantiene:
    segmentos_rep <-cbind(segmentos_rep, mcoin_object$mcoa$Tl1[,compX],mcoin_object$mcoa$Tl1[,compY]) 
    
    # Usaremos los elementos anteriores para dibujar el Scatter plot
    mcoa_ScorePlot <-plot_ly(as.data.frame(mcoin_object$mcoa$Tl1), 
                             x = ~mcoin_object$mcoa$Tl1[,compX], 
                             y = ~mcoin_object$mcoa$Tl1[,compY],
                             type ="scatter",
                             #hoverinfo = "text",  #Si la descomento para quitar coordenadas, las dos lineas de texto salen superpuestas. Se arregla quitando el /br y haciendo una linea larga, pero queda mas feo
                             text = paste("Sample:", do.call(rbind,lapply(strsplit(rownames(mcoin_object$mcoa$Tl1),"\\."), function(x) x[[1]])),
                                          "</br> OmicLevel: ",rep(datasetnames, each=nrow(mcoin_object$mcoa$SynVar))),
                             mode = "markers",
                             showlegend = TRUE,
                             color = factor(tratamientocolores,labels = c(unique(tratamientocolores))),
                             colors = paleta_tratamientos[1:length(unique(tratamientocolores))],
                             symbol = factor(rep(c(1:nrow(mcoin_object$mcoa$lambda)),each=nrow(mcoin_object$mcoa$SynVar)), labels=datasetnames),
                             symbols = simbolos[1:length(unique(tratamientocolores))], 
                             marker = list(size = 11)) #%>%
    mcoa_ScorePlot <- add_segments(mcoa_ScorePlot,
                                   x = ~segmentos_rep[,1], xend = ~segmentos_rep[,3],
                                   y = ~segmentos_rep[,2], yend = ~segmentos_rep[,4],
                                   name = NULL, 
                                   #mode = "line", 
                                   inherit = F, ## as? no hereda las movidas del anterior
                                   type = "lines", ## esta l?nea es la clave en verdad, el resto adorna.
                                   hoverinfo = "none",
                                   alpha = 0.75, 
                                   size = I(1),
                                   showlegend = FALSE,
                                   color = factor(tratamientocolores, #Luis: Añadi esto para colorear las lineas
                                                  labels = c(unique(tratamientocolores))),
                                   line = list(dash = "1px",  width = 2), ## TODO cambiar el color para que vaya en funci?n del nivel
                                   opacity = 1) %>%
      layout(xaxis = list(title = paste("Component",compX,sep=" "),titlefont=list(size=fontsizes[1]),tickfont=list(size=fontsizes[2])),
             yaxis = list(title = paste("Component",compY,sep=" "),titlefont=list(size=fontsizes[1]),tickfont=list(size=fontsizes[2])),
             title = "MCIA Co-Inertia Normed Scores Plot",
             titlefont = list(size=fontsizes[3]),
             legend = list(font = list(size=fontsizes[4])))
    return(mcoa_ScorePlot)
  }
  
  if(plottype=="Scree"){
    #Representamos los eigenvalues, y luego como sube la proporcion acumulada.
    
    # ***** Leyenda en ejes o en serires o ambas???? *********
    norm_eigenvalues <- mcoin_object$mcoa$pseudoeig/sum(mcoin_object$mcoa$pseudoeig) #Sacamos los autovalores, y os dividimos por la suma para tener el porcentaje de varianza en teoria recogida por cada tratamiento
    screeplotdata <- data.frame(Cs=paste("Comp",c(1:length(norm_eigenvalues)),sep=" "), ExplVar=norm_eigenvalues, Cumulat=cumsum(norm_eigenvalues)*100/sum(norm_eigenvalues),stringsAsFactors = FALSE) #Creamos tabla de datos con % varianza y % var acumulado
    screeplotdata$Cs <- factor(screeplotdata$Cs, levels = unique(screeplotdata$Cs)[order(screeplotdata$ExplVar, decreasing = TRUE)]) # Ordenamos los PCs segun valor
    
    mcoa_ScorePlot <- plot_ly(screeplotdata) %>%
      add_trace(x=~Cs,y=~ExplVar,type="bar",name="Eigenvalues",
                marker = list(color = '#C9EFF9')) %>%
      add_trace(x=~Cs, y=~Cumulat,type="scatter",mode="lines+markers",yaxis = "y2",name="Cumulative Eigenvalue",
                line = list(color = '#45171D')) %>%
      layout(title = "MCOA Screeplot", 
             xaxis = list(title = "",titlefont=list(size=fontsizes[1]),tickfont=list(size=fontsizes[2])),
             yaxis = list(side = 'left', title = "Eigenvalue", showgrid = F, zeroline = TRUE,titlefont=list(size=fontsizes[1]),tickfont=list(size=fontsizes[2])),
             yaxis2 = list(side = 'right', overlaying = "y", title = "Cumulative (Percent)", showgrid = FALSE, zeroline = FALSE,titlefont=list(size=fontsizes[1]),tickfont=list(size=fontsizes[2])),
             legend = list(x = 0.7, y = 0.1),
             titlefont = list(size=fontsizes[3]),
             legend = list(font = list(size=fontsizes[4])))
    
    return(mcoa_ScorePlot)
  }
  
  if(plottype=="Var"){
    if(levelchoice=="all"){
      if(useannot==FALSE){
    mcoa_ScorePlot <- plot_ly(as.data.frame(mcoin_object$mcoa$Tco), x = ~mcoin_object$mcoa$Tco[,compX], y = ~mcoin_object$mcoa$Tco[,compY],
                              text = paste("ID:", rownames(mcoin_object$mcoa$Tco),"</br> Description: ",annotations[,2],"</br> Mercator Bin: ",annotations[,4]), 
                              type="scatter",
                              mode = "markers", 
                              marker = list(size = 6, opacity=1), 
                              color=factor(rep(1:length(datasetnames),variablesxnivel), labels=datasetnames), 
                              colors =paleta_niveles[1:length(datasetnames)])
    
    mcoa_ScorePlot <- layout(mcoa_ScorePlot, title = "MCOA Variable Plot, all levels", 
                             xaxis = list(title = paste("Comp", compX, sep=" "),titlefont=list(size=fontsizes[1]),tickfont=list(size=fontsizes[2])),
                             yaxis = list(title = paste("Comp", compY, sep=" "),titlefont=list(size=fontsizes[1]),tickfont=list(size=fontsizes[2])),                                 titlefont = list(size=fontsizes[3]),
                             legend = list(font = list(size=fontsizes[4])))
    return(mcoa_ScorePlot)
    }
      if(useannot==TRUE){
      mapmanlevels <- as.numeric(unique(annotations[order(as.numeric(annotations[,3])),3]))
      mapmannames <- unique(annotations[order(as.numeric(annotations[,3])),4])
    mcoa_ScorePlot <- plot_ly(as.data.frame(mcoin_object$mcoa$Tco), 
                              x = ~mcoin_object$mcoa$Tco[,compX], 
                              y = ~mcoin_object$mcoa$Tco[,compY],
                              text = paste("ID:", rownames(mcoin_object$mcoa$Tco),"</br> Description: ",annotations[,2],"</br> Mercator Bin: ",annotations[,4]), 
                              type="scatter",
                              mode = "markers", 
                              marker = list(size = 6, opacity=1),
                              color = factor(annotations[,4], levels=mapmannames),
                              colors = paleta_mapman[mapmanlevels],
                              symbol =factor(rep(simbolos[1:length(datasetnames)],as.vector(mcoin_object$mcoa$blo)), labels=datasetnames), 
                              symbols = simbolos[1:length(datasetnames)])
    
    mcoa_ScorePlot <- layout(mcoa_ScorePlot, title = "MCOA Variable Plot, all levels", 
                             xaxis = list(title = paste("Comp", compX, sep=" "),titlefont=list(size=fontsizes[1]),tickfont=list(size=fontsizes[2])),
                             yaxis = list(title = paste("Comp", compY, sep=" "),titlefont=list(size=fontsizes[1]),tickfont=list(size=fontsizes[2])),
                             titlefont = list(size=fontsizes[3]),
                             legend = list(font = list(size=fontsizes[4])))
    return(mcoa_ScorePlot)
    }
  }
    if(levelchoice!="all"){
      if(levelchoice %in% datasetnames==FALSE)
        stop("Plot Var error: Please select a valid datasetname")
      # Definimos las filas a seleccionar para mostrar las variables correspondientes a cada dataset
      if(which(datasetnames==levelchoice)==1)
        filas <-c(1,mcoin_object$mcoa$blo[1])
      if(which(datasetnames==levelchoice)!=1)
        filas <- c((sum(mcoin_object$mcoa$blo[1:(which(datasetnames==levelchoice)-1)])+1), sum(mcoin_object$mcoa$blo[1:which(datasetnames==levelchoice)]))
      
      if(useannot==FALSE){
        mcoa_ScorePlot <- plot_ly(as.data.frame(mcoin_object$mcoa$Tco), 
                                  x = ~mcoin_object$mcoa$Tco[filas[1]:filas[2],compX], 
                                  y = ~mcoin_object$mcoa$Tco[filas[1]:filas[2],compY],
                                  text = paste("ID:", rownames(mcoin_object$mcoa$Tco)[filas[1]:filas[2]],"</br> Description: ",annotations[,2][filas[1]:filas[2]],"</br> Mercator Bin: ",annotations[,4][filas[1]:filas[2]]), 
                                  type="scatter",
                                  mode = "markers", 
                                  marker = list(size = 6, opacity=1))#antes estaba 0.85, 
        
        mcoa_ScorePlot <- layout(mcoa_ScorePlot, title = paste("MCOA Variable Plot, ",levelchoice," dataset",sep=""), 
                                 xaxis = list(title = paste("Comp", compX, sep=" "),titlefont=list(size=fontsizes[1]),tickfont=list(size=fontsizes[2])),
                                 yaxis = list(title = paste("Comp", compY, sep=" "),titlefont=list(size=fontsizes[1]),tickfont=list(size=fontsizes[2])),
                                 titlefont = list(size=fontsizes[3]),
                                 legend = list(font = list(size=fontsizes[4])))
        return(mcoa_ScorePlot)
      }
      if(useannot==TRUE){
        mapmanlevels <- as.numeric(unique(annotations[order(as.numeric(annotations[filas[1]:filas[2],3])),3]))
        mapmannames <- unique(annotations[order(as.numeric(annotations[filas[1]:filas[2],3])),4])
        mcoa_ScorePlot <- plot_ly(as.data.frame(mcoin_object$mcoa$Tco), 
                                  x = ~mcoin_object$mcoa$Tco[filas[1]:filas[2],compX], 
                                  y = ~mcoin_object$mcoa$Tco[filas[1]:filas[2],compY],
                                  text = paste("ID:", rownames(mcoin_object$mcoa$Tco)[filas[1]:filas[2]],"</br> Description: ",annotations[,2][filas[1]:filas[2]],"</br> Mercator Bin: ",annotations[,4][filas[1]:filas[2]]), 
                                  type="scatter",
                                  mode = "markers", 
                                  marker = list(size = 6, opacity=1),
                                  color = factor(annotations[filas[1]:filas[2],4], levels=mapmannames),
                                  colors = paleta_mapman[mapmanlevels])
    
        mcoa_ScorePlot <- layout(mcoa_ScorePlot, title = paste("MCOA Variable Plot, ",levelchoice," dataset",sep=""), 
                                 xaxis = list(title = paste("Comp", compX, sep=" "),titlefont=list(size=fontsizes[1]),tickfont=list(size=fontsizes[2])),
                                 yaxis = list(title = paste("Comp", compY, sep=" "),titlefont=list(size=fontsizes[1]),tickfont=list(size=fontsizes[2])),
                                 titlefont = list(size=fontsizes[3]),
                                 legend = list(font = list(size=fontsizes[4])))
        return(mcoa_ScorePlot)
      }
      
     }
  }
 
  if(plottype=="Biplot"){ 
    if(treatment!=4){  #Asi seleccionamos todas las visualizaciones que solo tengan un factor a mostrar.
      if(treatment==1) tratamientocolores <-treatmentstable[,1]
      if(treatment==2) tratamientocolores <-treatmentstable[,2]
      if(treatment==3) tratamientocolores <-treatmentstable[,3]
        
        if(levelchoice=="all"){
          if(useannot==FALSE){
            paletaintegrada <- c("#466fc7", paleta_tratamientos[1:length(unique(tratamientocolores))])
            mcoa_ScorePlot <- plot_ly() %>%
              add_trace(as.data.frame(mcoin_object$mcoa$Tco), 
                        x = ~mcoin_object$mcoa$Tco[,compX], 
                        y = ~mcoin_object$mcoa$Tco[,compY],
                        text = paste("ID:", rownames(mcoin_object$mcoa$Tco),"</br> Description: ",annotations[,2],"</br> Mercator Bin: ",annotations[,4]), 
                        type="scatter",
                        mode = "markers", 
                        marker = list(size = 5, opacity=1), 
                        color = as.factor("Loadings"),
                        colors = paletaintegrada,
                        symbol =factor(rep(simbolos[1:length(datasetnames)],as.vector(mcoin_object$mcoa$blo)), labels=datasetnames), 
                        symbols = simbolos[1:length(datasetnames)]) %>%
              add_trace(as.data.frame(mcoin_object$mcoa$SynVar), 
                      x = ~mcoin_object$mcoa$SynVar[,compX], 
                      y = ~mcoin_object$mcoa$SynVar[,compY],
                      text = paste("Sample:", rownames(mcoin_object$mcoa$SynVar), sep=" "),
                      showlegend = TRUE,
                      type="scatter",
                      mode = "markers", 
                      color = factor(tratamientocolores,labels = c(unique(tratamientocolores))),
                      colors = paletaintegrada[2:length(paletaintegrada)],
                      marker = list(size = 11),
                      yaxis="y",
                      xaxis="x") %>%

              layout(xaxis = list(title = paste("Component",compX,sep=" "),titlefont=list(size=fontsizes[1]),tickfont=list(size=fontsizes[2])),
                     yaxis = list(title = paste("Component",compY,sep=" "),titlefont=list(size=fontsizes[1]),tickfont=list(size=fontsizes[2])),
                     title = "MCIA Biplot",
                     titlefont = list(size=fontsizes[3]),
                     legend = list(font = list(size=fontsizes[4])))
       
            return(mcoa_ScorePlot)
          }
          if(useannot==TRUE){
            mapmanlevels <- as.numeric(unique(annotations[order(as.numeric(annotations[,3])),3]))
            mapmannames <- annotations[order(as.numeric(annotations[,3])),]
            mapmannames <- as.vector(unique(mapmannames[,4]))
            paletaintegrada <- c(paleta_mapman[mapmanlevels], paleta_tratamientos[1:length(unique(tratamientocolores))])
            

            mcoa_ScorePlot <- plot_ly() %>%
              
              add_trace(as.data.frame(mcoin_object$mcoa$Tco), 
                        x = ~mcoin_object$mcoa$Tco[,compX], 
                        y = ~mcoin_object$mcoa$Tco[,compY],
                        text = paste("ID:", rownames(mcoin_object$mcoa$Tco),"</br> Description: ",annotations[,2],"</br> Mercator Bin: ",annotations[,4]), 
                        type="scatter",
                        mode = "markers", 
                        marker = list(size = 6, opacity=1), 
                        color = factor(annotations[,4], levels=mapmannames),
                        colors = paletaintegrada,
                        symbol =factor(rep(simbolos[1:length(datasetnames)],as.vector(mcoin_object$mcoa$blo)), labels=datasetnames), 
                        symbols = simbolos[1:length(datasetnames)]) %>%
             
               add_trace(as.data.frame(mcoin_object$mcoa$SynVar), 
                        x = ~mcoin_object$mcoa$SynVar[,compX], 
                        y = ~mcoin_object$mcoa$SynVar[,compY],
                        # hoverinfo ="text",
                        text = paste("Sample:", rownames(mcoin_object$mcoa$SynVar), sep=" "),
                        showlegend = TRUE,
                        type="scatter",
                        mode = "markers", 
                        color = factor(tratamientocolores,labels = c(unique(tratamientocolores))),
                        colors = paleta_tratamientos[1:length(unique(tratamientocolores))],
                        marker = list(size = 11),
                        yaxis="y",
                        xaxis="x") %>%
               layout(xaxis = list(title = paste("Component",compX,sep=" "),titlefont=list(size=fontsizes[1]),tickfont=list(size=fontsizes[2])),
                     yaxis = list(title = paste("Component",compY,sep=" "),titlefont=list(size=fontsizes[1]),tickfont=list(size=fontsizes[2])),
                     title = "MCIA Biplot",
                     titlefont = list(size=fontsizes[3]),
                     legend = list(font = list(size=fontsizes[4])))
             return(mcoa_ScorePlot)
          }
        }
     
        if(levelchoice!="all"){ ## Revisar bien NO FUNCIONA
        if(levelchoice %in% datasetnames==FALSE)
          stop("Plot Var error: Please select a valid datasetname")
        # Definimos las filas a seleccionar para mostrar las variables correspondientes a cada dataset
        if(which(datasetnames==levelchoice)==1)
          filas <-c(1,mcoin_object$mcoa$blo[1])
        if(which(datasetnames==levelchoice)!=1)
          filas <- c((sum(mcoin_object$mcoa$blo[1:(which(datasetnames==levelchoice)-1)])+1), sum(mcoin_object$mcoa$blo[1:which(datasetnames==levelchoice)]))
        if(useannot==FALSE){
          paletaintegrada <- c("#466fc7", paleta_tratamientos[1:length(unique(tratamientocolores))])
          mcoa_ScorePlot <- plot_ly() %>%
            add_trace(as.data.frame(mcoin_object$mcoa$Tco), 
                      x = ~mcoin_object$mcoa$Tco[filas[1]:filas[2],compX], 
                      y = ~mcoin_object$mcoa$Tco[filas[1]:filas[2],compY],
                      text = paste("ID:", rownames(mcoin_object$mcoa$Tco)[filas[1]:filas[2]],"</br> Description: ",annotations[,2][filas[1]:filas[2]],"</br> Mercator Bin: ",annotations[,4][filas[1]:filas[2]]),
                      type="scatter",
                      mode = "markers", 
                      marker = list(size = 5, opacity=1), 
                      colors = paletaintegrada,
                      symbol =factor(rep(simbolos[1:length(datasetnames)],as.vector(mcoin_object$mcoa$blo)), labels=datasetnames), 
                      symbols = simbolos[1:length(datasetnames)]) %>%
            add_trace(as.data.frame(mcoin_object$mcoa$SynVar), 
                      x = ~mcoin_object$mcoa$SynVar[,compX], 
                      y = ~mcoin_object$mcoa$SynVar[,compY],
                      text = paste("Sample:", rownames(mcoin_object$mcoa$SynVar), sep=" "),
                      showlegend = TRUE,
                      type="scatter",
                      mode = "markers", 
                      color = factor(treatmentstable[,1],labels = c(unique(treatmentstable[,1]))),
                      colors = paletaintegrada[2:length(paletaintegrada)],
                      marker = list(size = 11),
                      yaxis="y",
                      xaxis="x") %>%

            layout(xaxis = list(title = paste("Component",compX,sep=" "),titlefont=list(size=fontsizes[1]),tickfont=list(size=fontsizes[2])),
                   yaxis = list(title = paste("Component",compY,sep=" "),titlefont=list(size=fontsizes[1]),tickfont=list(size=fontsizes[2])),
                   title = "MCIA Biplot",
                   titlefont = list(size=fontsizes[3]),
                   legend = list(font = list(size=fontsizes[4])))
            return(mcoa_ScorePlot)
          }
          if(useannot==TRUE){
            mapmanlevels <- as.numeric(unique(annotations[order(as.numeric(annotations[filas[1]:filas[2],3])),3]))
            mapmannames <- unique(annotations[order(as.numeric(annotations[filas[1]:filas[2],3])),4])
            mcoa_ScorePlot <- plot_ly() %>%
              add_trace(as.data.frame(mcoin_object$mcoa$SynVar), 
                        x = ~mcoin_object$mcoa$SynVar[,compX], 
                        y = ~mcoin_object$mcoa$SynVar[,compY],
                        # hoverinfo ="text",
                        text = paste("Sample:", rownames(mcoin_object$mcoa$SynVar), sep=" "),
                        showlegend = TRUE,
                        type="scatter",
                        mode = "markers", 
                        color = factor(tratamientocolores,labels = c(unique(tratamientocolores))),
                        colors = paleta_tratamientos[1:length(unique(tratamientocolores))],
                        #symbol = factor(rep(c(1:nrow(mcoin_object$mcoa$lambda)),each=nrow(mcoin_object$mcoa$SynVar)), labels=datasetnames),
                        symbols = simbolos[1:length(unique(tratamientocolores))], 
                        marker = list(size = 11),
                        yaxis="y",
                        xaxis="x") %>%
              add_trace(as.data.frame(mcoin_object$mcoa$Tco), 
                        x = ~mcoin_object$mcoa$Tco[filas[1]:filas[2],compX], 
                        y = ~mcoin_object$mcoa$Tco[filas[1]:filas[2],compY],
                        text = paste("ID:", rownames(mcoin_object$mcoa$Tco)[filas[1]:filas[2]],"</br> Description: ",annotations[,2][filas[1]:filas[2]],"</br> Mercator Bin: ",annotations[,4][filas[1]:filas[2]]), 
                        type="scatter",
                        mode = "markers", 
                        marker = list(size = 6, opacity=1),
                        color = factor(annotations[filas[1]:filas[2],4], levels=mapmannames),
                        colors = paleta_mapman[mapmanlevels]) %>%
              layout(xaxis = list(title = paste("Component",compX,sep=" "),titlefont=list(size=fontsizes[1]),tickfont=list(size=fontsizes[2])),
                     yaxis = list(title = paste("Component",compY,sep=" "),titlefont=list(size=fontsizes[1]),tickfont=list(size=fontsizes[2])),
                     title = "MCIA Biplot",
                     titlefont = list(size=fontsizes[3]),
                     legend = list(font = list(size=fontsizes[4])))
            return(mcoa_ScorePlot)
          }
        
        }
    }
    
    if(treatment==4){  #Asi seleccionamos las visualizaciones con 2 factores a mostrar
        if(levelchoice=="all"){
          if(useannot==FALSE){
            tratamientocolores <-treatmentstable[,1]
                      paletaintegrada <- c("#466fc7", paleta_tratamientos[1:length(unique(treatmentstable[,1]))])
  
            mcoa_ScorePlot <- plot_ly() %>%
              add_trace(as.data.frame(mcoin_object$mcoa$Tco), 
                        x = ~mcoin_object$mcoa$Tco[,compX], 
                        y = ~mcoin_object$mcoa$Tco[,compY],
                        text = paste("ID:", rownames(mcoin_object$mcoa$Tco),"</br> Description: ",annotations[,2],"</br> Mercator Bin: ",annotations[,4]), 
                        type="scatter",
                        mode = "markers", 
                        marker = list(size = 5, opacity=1), 
                        color = as.factor("Loadings"),
                        colors = paletaintegrada,
                        symbol =factor(rep(simbolos[1:length(datasetnames)],as.vector(mcoin_object$mcoa$blo)), labels=datasetnames), 
                        symbols = simbolos[(1+length(unique(treatmentstable[,2]))):(1+length(unique(treatmentstable[,2]))+length(unique(treatmentstable[,2])))]) %>%
                        
              add_trace(as.data.frame(mcoin_object$mcoa$SynVar), 
                        x = ~mcoin_object$mcoa$SynVar[,compX], 
                        y = ~mcoin_object$mcoa$SynVar[,compY],
                        text = paste("Sample:", rownames(mcoin_object$mcoa$SynVar), sep=" "),
                        showlegend = TRUE,
                        type="scatter",
                        mode = "markers", 
                        color = factor(treatmentstable[,1],labels = c(unique(treatmentstable[,1]))),
                        colors = paletaintegrada[2:length(paletaintegrada)],
                        marker = list(size = 11),
                        symbol = factor(treatmentstable[,2], labels=unique(treatmentstable[,2])),
                        symbols = simbolos[1:length(unique(treatmentstable[,2]))], 
                        yaxis="y",
                        xaxis="x") %>%
              
              layout(xaxis = list(title = paste("Component",compX,sep=" "),titlefont=list(size=fontsizes[1]),tickfont=list(size=fontsizes[2])),
                     yaxis = list(title = paste("Component",compY,sep=" "),titlefont=list(size=fontsizes[1]),tickfont=list(size=fontsizes[2])),
                     title = "MCIA Biplot",
                     titlefont = list(size=fontsizes[3]),
                     legend = list(font = list(size=fontsizes[4])))
            return(mcoa_ScorePlot)
          }

          if(useannot==TRUE){
            mapmanlevels <- as.numeric(unique(annotations[order(as.numeric(annotations[,3])),3]))
            mapmannames <- annotations[order(as.numeric(annotations[,3])),]
            mapmannames <- as.vector(unique(mapmannames[,4]))
            paletaintegrada <- c(paleta_mapman[mapmanlevels], paleta_tratamientos[1:length(unique(treatmentstable[,1]))])
            
            mcoa_ScorePlot <- plot_ly() %>%
              add_trace(as.data.frame(mcoin_object$mcoa$Tco), 
                        x = ~mcoin_object$mcoa$Tco[,compX], 
                        y = ~mcoin_object$mcoa$Tco[,compY],
                        text = paste("ID:", rownames(mcoin_object$mcoa$Tco),"</br> Description: ",annotations[,2],"</br> Mercator Bin: ",annotations[,4]), 
                        type="scatter",
                        mode = "markers", 
                        marker = list(size = 5, opacity=1), 
                        color = as.factor("Loadings"),
                        colors = paletaintegrada,
                        symbol =factor(rep(simbolos[1:length(datasetnames)],as.vector(mcoin_object$mcoa$blo)), labels=datasetnames), 
                        symbols = simbolos[1:length(datasetnames)]) %>%
              
              add_trace(as.data.frame(mcoin_object$mcoa$SynVar), 
                        x = ~mcoin_object$mcoa$SynVar[,compX], 
                        y = ~mcoin_object$mcoa$SynVar[,compY],
                        text = paste("Sample:", rownames(mcoin_object$mcoa$SynVar), sep=" "),
                        showlegend = TRUE,
                        type="scatter",
                        mode = "markers", 
                        color = factor(treatmentstable[,1],labels = c(unique(treatmentstable[,1]))),
                        colors = paletaintegrada[2:length(paletaintegrada)],
                        marker = list(size = 11),
                        symbol = factor(treatmentstable[,2], labels=unique(treatmentstable[,2])),
                        symbols = simbolos[1:length(unique(treatmentstable[,2]))], 
                        yaxis="y",
                        xaxis="x") %>%

              layout(xaxis = list(title = paste("Component",compX,sep=" "),titlefont=list(size=fontsizes[1]),tickfont=list(size=fontsizes[2])),
                     yaxis = list(title = paste("Component",compY,sep=" "),titlefont=list(size=fontsizes[1]),tickfont=list(size=fontsizes[2])),
                     title = "MCIA Biplot",
                     titlefont = list(size=fontsizes[3]),
                     legend = list(font = list(size=fontsizes[4])))
            return(mcoa_ScorePlot)
          }
        }
        
      
        if(levelchoice!="all"){
          if(levelchoice %in% datasetnames==FALSE)
            stop("Plot Var error: Please select a valid datasetname")
          # Definimos las filas a seleccionar para mostrar las variables correspondientes a cada dataset
          if(which(datasetnames==levelchoice)==1)
            filas <-c(1,mcoin_object$mcoa$blo[1])
          if(which(datasetnames==levelchoice)!=1)
            filas <- c((sum(mcoin_object$mcoa$blo[1:(which(datasetnames==levelchoice)-1)])+1), sum(mcoin_object$mcoa$blo[1:which(datasetnames==levelchoice)]))
          
          if(useannot==FALSE){
            mcoa_ScorePlot <- plot_ly() %>%
              add_trace(as.data.frame(mcoin_object$mcoa$SynVar), 
                        x = ~mcoin_object$mcoa$SynVar[,compX], 
                        y = ~mcoin_object$mcoa$SynVar[,compY],
                        # hoverinfo ="text",
                        text = paste("Sample:", rownames(mcoin_object$mcoa$SynVar), sep=" "),
                        showlegend = TRUE,
                        type="scatter",
                        mode = "markers", 
                        color = factor(treatmentstable[,1], labels = unique(treatmentstable[,1])),
                        colors = paleta_tratamientos[1:length(unique(treatmentstable[,1]))],
                        symbol = factor(treatmentstable[,2], labels=unique(treatmentstable[,2])),
                        symbols = simbolos[1:length(unique(treatmentstable[,2]))], 
                        marker = list(size = 11),
                        yaxis="y",
                        xaxis="x") %>%
              add_trace(as.data.frame(mcoin_object$mcoa$Tco), 
                        x = ~mcoin_object$mcoa$Tco[filas[1]:filas[2],compX], 
                        y = ~mcoin_object$mcoa$Tco[filas[1]:filas[2],compY],
                        text = paste("ID:", rownames(mcoin_object$mcoa$Tco)[filas[1]:filas[2]],"</br> Description: ",annotations[,2][filas[1]:filas[2]],"</br> Mercator Bin: ",annotations[,4][filas[1]:filas[2]]),
                        type="scatter",
                        mode = "markers", 
                        marker = list(size = 5, opacity=1), 
                        color=factor(rep(1:length(datasetnames),variablesxnivel), labels=datasetnames), 
                        colors =paleta_niveles[1:length(datasetnames)]) %>%
              layout(xaxis = list(title = paste("Component",compX,sep=" "),titlefont=list(size=fontsizes[1]),tickfont=list(size=fontsizes[2])),
                     yaxis = list(title = paste("Component",compY,sep=" "),titlefont=list(size=fontsizes[1]),tickfont=list(size=fontsizes[2])),
                     title = "MCIA Biplot",
                     titlefont = list(size=fontsizes[3]),
                     legend = list(font = list(size=fontsizes[4])))
            
            return(mcoa_ScorePlot)
          }
          if(useannot==TRUE){
            mapmanlevels <- as.numeric(unique(annotations[order(as.numeric(annotations[filas[1]:filas[2],3])),3]))
            mapmannames <- unique(annotations[order(as.numeric(annotations[filas[1]:filas[2],3])),4])
            mcoa_ScorePlot <- plot_ly() %>%
              add_trace(as.data.frame(mcoin_object$mcoa$SynVar), 
                        x = ~mcoin_object$mcoa$SynVar[,compX], 
                        y = ~mcoin_object$mcoa$SynVar[,compY],
                        # hoverinfo ="text",
                        text = paste("Sample:", rownames(mcoin_object$mcoa$SynVar), sep=" "),
                        showlegend = TRUE,
                        type="scatter",
                        mode = "markers", 
                        color = factor(treatmentstable[,1], labels = unique(treatmentstable[,1])),
                        colors = paleta_tratamientos[1:length(unique(treatmentstable[,1]))],
                        symbol = factor(treatmentstable[,2], labels=unique(treatmentstable[,2])),
                        symbols = simbolos[1:length(unique(treatmentstable[,2]))], 
                        marker = list(size = 11),
                        yaxis="y",
                        xaxis="x") %>%
              add_trace(as.data.frame(mcoin_object$mcoa$Tco), 
                        x = ~mcoin_object$mcoa$Tco[filas[1]:filas[2],compX], 
                        y = ~mcoin_object$mcoa$Tco[filas[1]:filas[2],compY],
                        text = paste("ID:", rownames(mcoin_object$mcoa$Tco)[filas[1]:filas[2]],"</br> Description: ",annotations[,2][filas[1]:filas[2]],"</br> Mercator Bin: ",annotations[,4][filas[1]:filas[2]]), 
                        type="scatter",
                        mode = "markers", 
                        marker = list(size = 6, opacity=1),
                        color = factor(annotations[filas[1]:filas[2],4], levels=mapmannames),
                        colors = paleta_mapman[mapmanlevels]) %>%
              layout(xaxis = list(title = paste("Component",compX,sep=" "),titlefont=list(size=fontsizes[1]),tickfont=list(size=fontsizes[2])),
                     yaxis = list(title = paste("Component",compY,sep=" "),titlefont=list(size=fontsizes[1]),tickfont=list(size=fontsizes[2])),
                     title = "MCIA Biplot",
                     titlefont = list(size=fontsizes[3]),
                     legend = list(font = list(size=fontsizes[4])))
            return(mcoa_ScorePlot)
          }
          
          
        }
        
      }
       
    }
        
  if(plottype=="Pseudoeigen"){
    #Representamos como de proximos o alejados están los niveles por su autovalor.
    
    mcoa_ScorePlot <- plot_ly(as.data.frame(mcoin_object$mcoa$cov2), x = ~mcoin_object$mcoa$cov2[,compX], y = ~mcoin_object$mcoa$cov2[,compY],
                              type="scatter",
                              mode = "markers", 
                              marker = list(size = 11, opacity=1), 
                              color=factor(c(1:length(datasetnames)), labels=datasetnames), 
                              colors =paleta_niveles[1:length(datasetnames)], 
                              hoverinfo="none")
    
    mcoa_ScorePlot <- layout(mcoa_ScorePlot, title = "Pseudoeigenvalues. All datasets", 
                             xaxis = list(title = paste("Pseudoeigenvalue", compX, sep=" "),titlefont=list(size=fontsizes[1]),tickfont=list(size=fontsizes[2])),
                             yaxis = list(title = paste("Pseudoeigenvalue", compY, sep=" "),titlefont=list(size=fontsizes[1]),tickfont=list(size=fontsizes[2])),
                             titlefont = list(size=fontsizes[3]),
                             legend = list(font = list(size=fontsizes[4])))
    return(mcoa_ScorePlot)
  }
  
  if(plottype=="Composite"){
    if(!is.null(compZ)) stop("Composite plot can only draw 2D plots. compZ must be null") 
      #Representamos los eigenvalues, y luego como sube la proporcion acumulada.
    p1 <- mcia_plot(mcoin_list,treatment, compX, compY, plottype = "SynVar")
    p2 <- mcia_plot(mcoin_list,treatment, compX, compY, plottype = "Var", useannot=FALSE,levelchoice="all")
    p3 <- mcia_plot(mcoin_list,treatment, compX, compY, plottype = "CoIn")
    p4 <- mcia_plot(mcoin_list,treatment, compX, compY, plottype = "Pseudoeigen")
    mcoa_ScorePlot <- subplot(subplot(p1,p2),subplot(p3,p4),nrows=2,shareX = F,shareY = F,margin = 0.04)
    mcoa_ScorePlot <- layout(mcoa_ScorePlot, title = "MCOA Analysis")
        return(mcoa_ScorePlot)
  }
  
  if(plottype=="All"){
    #Aqui habra que crear una llamada recursiva a todos los graficos posibles y devolver una lista.
    #tambien recursivo pa los distintos niveles.
    if(!is.null(compZ)) stop("Please select only two dimensions for plotting all figures")
    cat("EXPORTING ANALYSIS PLOTS\nA single core will be used. It may take a while...\n ")
    cat("Generating and saving SynVar plot. Step 1/6\n")
    export_plot(mcia_plot(mcoin_list,treatment=treatment, compX=compX, compY=compY, plottype="SynVar",useannot=useannot,levelchoice=levelchoice),"mcoa_synvar.pdf")
    cat("Generating and saving Scree plot. Step 2/6\n")
    export_plot(mcia_plot(mcoin_list,treatment=treatment, compX=compX, compY=compY, plottype="Scree",useannot=useannot,levelchoice=levelchoice),"mcoa_scree.pdf")
    cat("Generating and saving Var plot. Step 3/6\n")
    export_plot(mcia_plot(mcoin_list,treatment=treatment, compX=compX, compY=compY, plottype="Var",useannot=useannot,levelchoice=levelchoice),"mcoa_var.pdf")
    cat("Generating and saving CoInertia plot. Step 4/6\n")
    export_plot(mcia_plot(mcoin_list,treatment=treatment, compX=compX, compY=compY, plottype="CoIn",useannot=useannot,levelchoice=levelchoice),"mcoa_coin.pdf")
    cat("Generating and saving Var plot. Step 5/6\n")
    export_plot(mcia_plot(mcoin_list,treatment=treatment, compX=compX, compY=compY, plottype="Pseudoeigen",useannot=useannot,levelchoice=levelchoice),"mcoa_pseudoeigen.pdf")
    cat("Generating and saving Composite plot. Step 6/6\n")
    export_plot(mcia_plot(mcoin_list,treatment=treatment, compX=compX, compY=compY, plottype="Composite",useannot=useannot,levelchoice=levelchoice),"mcoa_composite.pdf")
  }
}
