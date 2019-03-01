#k-means script almost finished

#Bueno esto va como siempre... initial row, initial column = coordenadas de donde empieza la parte numerica de la matriz;
#treatment column = columna que nos indica de que tratamiento es replica cada fila
#Original matrix structure: variables en columnas, replicas/muestras en filas
#scalation es un argumento cutre que he metido para escalar los valores o lo que sea, lo he puesto numerico
  #scalation=1 escalar y NO centrar 
  #scalation=2 escalar y centrar
  #scalation=3 escala entre 0 y 1 por columna
  #scalation=4 escala entre 0 y 1 por fila

kmeans_analysis<-function(listadatos,annotation=NULL,initialrow=1, initialcolumn=3, treatment1col=1, treatment2col=2, treatment=1,omiclevel="all",scalation,clusters=c(10,20),show.elbow.plot=TRUE){
  require(plyr)
  options(warn=-1)
  if(class(listadatos)!="POL") stop("A POL object is required. Please run pre-processing steps")
  datasetnames<-names(listadatos)
  
  #Generamos un vector para hacer el split por tratamientos
  if(treatment==1) vectortratamientos <- as.factor(listadatos[[1]][initialrow:nrow(listadatos[[1]]),treatment1col])
  if(treatment==2) vectortratamientos <- as.factor(listadatos[[1]][initialrow:nrow(listadatos[[1]]),treatment2col])
  if(treatment==3) vectortratamientos <- as.factor(paste(listadatos[[1]][initialrow:nrow(listadatos[[1]]),treatment1col],listadatos[[1]][initialrow:nrow(listadatos[[1]]),treatment2col], sep="&"))
  
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
  datos <- split(as.data.frame(listaintermedia),vectortratamientos,drop=T)
  datos<-datos[unique(vectortratamientos)]
  #a veces en R al editar o generar una matriz los valores numericos se convierten en factores, para volver a valores numericos
  #los defines primero como caracteres y luego como numero, si pones as.numeric  a pelo MAL MuY MAL, porque te devuelve el ordinal 
  #de cada factor, es decir te los ordena de menor a mayor y el numero que te devuelve es la posicion en la lista
  #no da error, y no te das cuenta...
  datos<-lapply(datos, function(datos) t(apply(datos,1,function(datos) as.numeric(as.character(datos)))))
  #media de las replicas por cada variable
  treat.means<-lapply(datos, colMeans)
  #juntarlo en una tabla, gracias al unique no se van a descolocar nombres
  means.matrix<-do.call("rbind",treat.means)
  #y les volvemos a poner los nombres a las variables 
  colnames(means.matrix)<-colnames(listaintermedia)
  #hasta aqui la matriz para el k-means clustering
  #dim(means.matrix)=n? treatments x n? variables
  #Different choices to scale/center/percentage... whatever
    #scale
  if(scalation==1){meansfmatrix<-scale(means.matrix,center=F)}
    #scale and center
  else if (scalation==2){meansfmatrix<-scale(means.matrix)}
    #percentaje (tanto por 1) por columna, a mi me gusta este!
  else if (scalation==3){meansfmatrix<-sweep(means.matrix,2,colSums(means.matrix),"/" )  }
    #percentaje (tanto por 1) por fila, otra opcion
  else if (scalation==4){meansfmatrix<-sweep(means.matrix,1,rowSums(means.matrix),"/" )  }
 
  ## Ahora que los datos estan bien nos ponemos con el kmeans
  grupos<-c()
  grupos.ensayo <-c()
  elbow.data<-c()
  lista.kmeans <-list()
  kmeans_matrix <- t(meansfmatrix)
  if(clusters=="auto"){} #sin implementar aqui. La idea es rango de kmeans entre 5 y 50, y seleccionar rango optimo
  # hacemos el analisis kmeans
  j=1
  for(i in clusters[1]:clusters[2]){
    kmeans.result <-kmeans_helper(kmeans_matrix, i)
    grupos.ensayo <- cbind(grupos.ensayo,kmeans.result[[1]])
    colnames(grupos.ensayo)[j] <- paste(i,"Clusters",sep=" ")
    elbow.data <- c(elbow.data,kmeans.result[[2]])
    names(elbow.data)[j] <-paste(i,"Clusters",sep=" ")
    lista.kmeans[[j]] <-as.data.frame(kmeans.result[[3]])
    names(lista.kmeans)[[j]] <- paste(i,"Clusters",sep=" ")
    j=j+1
  }
  if(show.elbow.plot==TRUE){
    plot(elbow.data, type="o", xlab="Number of clusters, K", ylab="Cumulative within-clusters sum of squares",xaxt = "n")
    axis(1, at=1:length(elbow.data), labels=names(elbow.data))
  } 
  resultado <-list(meansfmatrix,annotmatrix,lista.kmeans,grupos.ensayo,elbow.data)
  names(resultado) <- c("kmeans_matrix","annotations","kmeans_list","kmeans_group","elbow.data")
  class(resultado) <- "kmeansanalysis"
   return(resultado)
}

kmeans_helper<-function(matrix,i){
  set.seed(5881)
  kmeans.result <- kmeans(matrix, i,iter.max=30) # Hacemos un k.means clustering (see below)
  grupos <- kmeans.result$cluster #sacamos a que grupo pertenece cada variable
  withinss <- sum(kmeans.result$withinss)
  datosggplot <- as.data.frame(cbind(matrix, grupos))
  datosggplot <- cbind(datosggplot, rownames(matrix))
  colnames(datosggplot)[ncol(datosggplot)] <- "ID"
 result <- list(grupos,withinss,datosggplot)
 return(result)
} 
  
kmeans_plot <- function(kmeansobject,clusternumber=NULL,fontsizes=c(14,10,16,12)){
  require(plyr)
  require(reshape2)
  require(ggplot2)
  require(plotly)
  options(warn=-1)
  if(class(kmeansobject)!="kmeansanalysis") stop("A kmeansanalysis object is required. Please run kmeans_analysis  first")
  if(!is.null(clusternumber)){
    clusteringsdone<-as.numeric(gsub(" Clusters","",names(kmeansobject$kmeans_list),fixed = TRUE)) #sacamos los numeros de los clusters hechos
    if(clusternumber%in%clusteringsdone==FALSE) stop("No se hizo kmeans con tantos clusters. repetir") #vemos que ese numero de grupos se haya hecho de forma efectiva
    cluster<-which(as.numeric(gsub(" Clusters","",names(kmeansobject$kmeans_list),fixed = TRUE))==clusternumber) #sacamos que posicion en la lista de kmeans list ocupa el numero de clusters deseado
    datosggplot<-kmeansobject$kmeans_list[[cluster]] #sacamos los datos requeridos
    datosggplot <- melt(datosggplot, (ncol(datosggplot)-1):(ncol(datosggplot))) #melt para ponerlo en formato para ggplot2
    datosggplot[,1] <- paste("Cluster", datosggplot[,1], sep= " ") #cambiamos los nombres de gurpos para poner cluster antes y que salga bonito
    myplot <- ggplot(datosggplot, aes(x=variable, y=value, group=ID, colour=as.factor(grupos), text=datosggplot$ID)) +
      facet_wrap(~ grupos,scales="free") + 
      geom_line(linetype= "dotted",alpha = 1/6) +
      stat_summary(aes(group = grupos), fun.y = mean, geom = "line") + 
      theme_minimal()+ theme(legend.position="none",axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
    
    gg<-ggplotly(myplot, tooltip = c("text"))
    
    gg<-layout(gg, title = "K-mean grouping", titlefont = list(size=fontsizes[3]),legend = list(font = list(size=fontsizes[4])))
    #class(gg)<-c("plotly","htmlwidget","singlekmeansplot")
    
    singlekmeansplot <-list(gg,myplot)
    class(singlekmeansplot) <- "singlekmeansplot"
    return(singlekmeansplot)
    }
  #ahora aqui exportaremos la lista de graficos
  ggplotlist <- c()
  myplotlist <- c()
  cat("Generating k-means plots. It may take a while\n")
  pb <- txtProgressBar(min = 0, max = length(kmeansobject$kmeans_list), style = 3)
  for(i in 1:length(kmeansobject$kmeans_list)){
    datosggplot<-kmeansobject$kmeans_list[[i]] #sacamos los datos requeridos
    datosggplot <- melt(datosggplot, (ncol(datosggplot)-1):(ncol(datosggplot))) #melt para ponerlo en formato para ggplot2
    datosggplot[,1] <- paste("Cluster", datosggplot[,1], sep= " ") #cambiamos los nombres de gurpos para poner cluster antes y que salga bonito
    myplot <- ggplot(datosggplot, aes(x=variable, y=value, group=ID, colour=as.factor(grupos), text=datosggplot$ID)) +
      facet_wrap(~ grupos,scales="free") + 
      geom_line(linetype= "dotted",alpha = 1/6) +
      stat_summary(aes(group = grupos), fun.y = mean, geom = "line") + 
      theme_minimal()+ theme(legend.position="none",axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
    
    gg<-ggplotly(myplot, tooltip = c("text"))
    
    gg<-layout(gg, title = "K-mean grouping", titlefont = list(size=fontsizes[3]),legend = list(font = list(size=fontsizes[4])))
    myplotlist[[i]] <- myplot #esto es una GRAN chapuza, pero no se exportar de otra forma
    ggplotlist[[i]] <- gg
    setTxtProgressBar(pb, i)
   
  }
  names(ggplotlist)<-names(kmeansobject$kmeans_list)
  names(myplotlist)<-names(kmeansobject$kmeans_list)
  #class(ggplotlist)<-c("plotly","htmlwidget","multiplekmeansplot")
  kmeansmultiple<-list(ggplotlist,myplotlist)
  class(kmeansmultiple)<-"kmeansmultipleplot"
  return(kmeansmultiple)   
  }
