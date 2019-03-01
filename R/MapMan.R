# Plottype = pie no funciona.
# Error in add_trace(p, ...) : objeto 'mapman' no encontrado 
# 11. add_trace(p, ...) 
#10. add_trace_classed(p, class = "plotly_pie", values = values, labels = labels, 
#                  type = "pie", ..., data = data, inherit = inherit) 
#9.add_pie(., mapman, labels = ~rownames(mapman), values = ~mapman[, 
#                                                                1], name = colnames(mapman)[1], marker = list(colors = paleta_mapman), 
#        domain = list(x = c(0, 0.4), y = c(0.4, 1)), sort = FALSE) 
#8.function_list[[i]](value) 
#7.freduce(value, `_function_list`) 
#6.`_fseq`(`_lhs`) 
#5.eval(quote(`_fseq`(`_lhs`)), env, env) 
#4.eval(quote(`_fseq`(`_lhs`)), env, env) 
#3.withVisible(eval(quote(`_fseq`(`_lhs`)), env, env)) 
#2.plot_ly() %>% add_pie(mapman, labels = ~rownames(mapman), values = ~mapman[, 
#                                                                           1], name = colnames(mapman)[1], marker = list(colors = paleta_mapman), 
#                      domain = list(x = c(0, 0.4), y = c(0.4, 1)), sort = FALSE) %>% 
#  add_pie(mapman, labels = ~rownames(mapman), values = ~mapman[,  ... at f.R#650
#1. mapman_plot(gruposmapman_prot, plottype = "pie") 
#                                                               

mapman_group<-function(listadatos, annotation, initialrow=1, initialcolumn=3, treatment1col=1, treatment2col=2, treatment=1, method="sum", group=TRUE, omiclevel=NULL, threshold=c(0,0), stats=FALSE,scalevalues=TRUE,normbyvariables=FALSE){
  options(warn=-1)
  require(dplyr)

  # Checks iniciales
  if(class(listadatos)!="POL") stop("\nA POL object is required. Please run pre-processing steps")
  if(is.null(annotation)) stop("\nAnnotation matrix is required. Please run MapMan to annotate your variables\nVisit https://mapman.gabipd.org for more information.")
  datasetnames<-names(listadatos)
  #if(!is.null(omiclevel)) if(omiclevel %in% datasetnames==FALSE) stop("Error: Please select a valid datasetname")

  #Nos aseguramos de quitar las columnas vacias
  listadatos <- lapply(listadatos, function(x) RemoveEmptyColumns(x, initialrow=initialrow,initialcolumn=initialcolumn))

  #Extraemos los valores significativos de las tablas
  listadatos2 <- lapply(listadatos, function(x)   x[initialrow:nrow(x),initialcolumn:ncol(x)])


  ##
  ## Apano para crear una tabla de identificaciones a usar posteriormente. Solo MapMan.
  ## es un poco chapuzas pero bueno, funciona
  ##
  #
  # Importante: chequear que la tabla de identificaciones no tenga duplicados, sino es un lio porque el dplyr crea filas extra y no cuadra
  # para los dibujos. Pensar en algo para que si hay duplicados de error y no siga.

  # Seleccionamos los datos
  if(is.null(omiclevel)) {
    variablenames <- lapply(listadatos2, function(x) colnames(x))
    variablenames <-unlist(variablenames)
    listaintermedia <-do.call(cbind,listadatos2)  #### Esto es inconsistente. Cambia nombres de columnas, a veces, nombretabla.nombrevar, otras veces solo nombrevar
    colnames(listaintermedia) <- variablenames #### para arreglar lode arriba entre tanto.
  }
  if(!is.null(omiclevel)) {
    sapply(omiclevel, function(x) if(x %in% datasetnames==FALSE) stop("Error: Please select a valid datasetname"))

    if(length(omiclevel)>1){  listaintermedia <-do.call(cbind,listadatos2[omiclevel])
    variablenames <- lapply(listadatos2[omiclevel], function(x) colnames(x))
    variablenames <-unlist(variablenames)
    colnames(listaintermedia) <- variablenames}

    if(length(omiclevel)==1) listaintermedia <-listadatos2[[omiclevel]]
    variablenames <-colnames(listaintermedia)
  }

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

  annotmatrix <-as.matrix(annotmatrix)
  annotmatrix <- annotmatrix[order(as.numeric(annotmatrix[,3])),]
  grupos<-sapply(unique(annotmatrix[,4]), function(x) annotmatrix[annotmatrix[,4]==x,1])

  # Generamos una matriz binaria, para contar el numero de elementos por categoria. Definimos valor umbral para ser considerada o no
  matrix_number_elements <- listaintermedia
  for(i in 1:nrow(matrix_number_elements)){
    for(j in 1:ncol(matrix_number_elements)){
      if(matrix_number_elements[i,j]>=threshold[1]){
        matrix_number_elements[i,j]<-1
      }
      else {
        matrix_number_elements[i,j]<-0
      }
    }
  }
  matrix_number_elements<-as.data.frame(matrix_number_elements)

  # Generamos una matriz con las sumas de las variables correspondientes a cada categoria mapman
  listagrupada<- sumpergroup(listaintermedia,grupos)


  # Definimos el vector para dividir los datos en funcion del tratamiento
  if(treatment==1) splitvector <- as.vector(listadatos[[1]][,treatment1col])
  if(treatment==2) splitvector <- as.vector(listadatos[[1]][,treatment2col])
  if(treatment==3) splitvector <- as.vector(paste(listadatos[[1]][,treatment1col], listadatos[[1]][,treatment2col],sep="&"))

  # Definimos la lista con los distintos tratamientos.
  listagrupada<-as.data.frame(listagrupada)
  result<-listagrupada
  listaxtratamientos <- split(listagrupada,splitvector,drop=T)

  # Generamos una matriz con las sumas del numero de variables correspondientes a cada categoria mapman

  matrix_number_elements_treatments<- as.data.frame(sumpergroup(matrix_number_elements,grupos))


  # Si se agrupa por tratamiento
  if(group==TRUE){

    medias<-lapply(listaxtratamientos, colMeans)
    medias<-as.data.frame(do.call(cbind, medias))

    mediassd <- lapply(listaxtratamientos, function(x) apply(x,2,sd))
    mediassd<-as.data.frame(do.call(cbind, mediassd))

    #
    # NUMERO DE ELEMENTOS POR CATEGORIA Y TRATAMIENTO
    # Sacamos el numero de variables por categoria
    #

    elementos_categoriaxtratamientos <-split(matrix_number_elements,splitvector,drop=T)
    elementos_categoriaxtratamientos <-lapply(elementos_categoriaxtratamientos, colMeans)
    elementos_categoriaxtratamientos <-as.data.frame(do.call(cbind,elementos_categoriaxtratamientos))
    elementos_categoriaxtratamientos[elementos_categoriaxtratamientos<threshold[2]]<-0 #Aqui metemos criterio de consistencia
    matrix_number_elements_treatments<- as.data.frame(sumpergroup(t(elementos_categoriaxtratamientos),grupos))

    #
    # NUMERO DE ELEMENTOS POR CATEGORIA Y TRATAMIENTO
    # Sacamos el numero de variables por categoria
    #
    if(stats==TRUE){
      a= c()
      pvalues=c()
      if(normbyvariables==TRUE) {
      listagrupada <- listagrupada/t(matrix_number_elements_treatments)
      }
      if(scalevalues==TRUE) listagrupada <- scale(listagrupada)
      #hacemos anova
      datamatrix <- as.data.frame(listagrupada)
      i=1
      while(i<=ncol(datamatrix)){
        res.anova=aov(datamatrix[,i]~as.factor(splitvector),data=datamatrix,na.rm=TRUE)
        pvalues <- c(pvalues,anova(res.anova)$'Pr(>F)'[1:1])
        i=i+1
        names(pvalues)[1:length(pvalues)] <- colnames(datamatrix)
      }
      result <-list(medias,mediassd,t(matrix_number_elements_treatments),pvalues,splitvector)
      names(result) <- c("mean","sd","number_elements_category","pvalues","treatments")
      class(result) <- "MMOG"
      return(result)
    }

    result <-list(medias,mediassd,t(matrix_number_elements_treatments),splitvector)
    names(result) <- c("mean","sd","number_elements_category","treatments")
    class(result) <- "MMOG"
    return(result)
  }

  result <- list(t(listagrupada),t(matrix_number_elements_treatments),splitvector)
  names(result) <- c("mean","number_elements_category","treatments")
  class(result) <- "MMO"

  return(result)

}

sumpergroup <-function(lista,grupos){
  p1<-lapply(grupos,function(x) as.data.frame(lista[,colnames(lista)%in%x]))
  p2 <-lapply(p1, function(x) rowSums(x,na.rm=T))
  p3 <- do.call(cbind,p2)
  return(p3)
}

mapman_plot <- function(listaMMO, plottype="bar", ploterror=TRUE, logscale=FALSE, hmparameters=c("manhattan","manhattan","ward.D",FALSE,"row"),normbyelements=FALSE,fontsizes=c(14,10,16,12)){
  require(plotly)
  if(class(listaMMO) %in% c("MMOG","MMO")==FALSE) stop("Not valid object. Please run mapman_group first")
  values <- listaMMO$mean
  sd <- listaMMO$sd
  treatments <- factor(listaMMO$treatments,labels=unique(listaMMO$treatments))
  if(class(listaMMO)=="MMO") ploterror=FALSE
  if(normbyelements==TRUE){
    if(class(listaMMO)=="MMOG") values <-listaMMO$mean/listaMMO$number_elements_category
  }

  if(plottype=="bar"){
    if(ploterror==T){
      p<-plot_ly(as.data.frame(values),
                 x=rownames(values),
                 y = values[,1],
                 type = "bar",
                 color = colnames(values)[1],
                 colors = paleta_tratamientos[1:ncol(values)],
                 showlegend = TRUE,
                 error_y=list(value= sd[1],color = "#656a73",thickness=1, width=1, symmetric=TRUE))
      if(ncol(values)>1){
        for(i in 2:ncol(values)){
          p<-add_trace(p,
                       x=rownames(values),
                       y = values[,i],
                       type = "bar",
                       color=colnames(values)[i],
                       error_y=list(value= sd[i],color = "#656a73",thickness=1, symmetric=TRUE))
        }
      }
    }

    if(ploterror==FALSE){
      p<-plot_ly(as.data.frame(values),
                 x=rownames(values),
                 y = values[,1],
                 type = "bar",
                 color = colnames(values)[1],
                 colors = paleta_tratamientos[1:ncol(values)],
                 showlegend = TRUE)
      if(ncol(values)>1){
        for(i in 2:ncol(values)){
          p<-add_trace(p,
                       x=rownames(values),
                       y = values[,i],
                       type = "bar",
                       color=colnames(values)[i])


        }
      }
    }

    p<-layout(p,
              margin = list(b = 160),
              xaxis=list(tickangle = 45,titlefont=list(size=fontsizes[1]),tickfont=list(size=fontsizes[2])),
              yaxis = list(title = "",type = "",tickformat="e",titlefont=list(size=fontsizes[1]),tickfont=list(size=fontsizes[2])),
              title = "Abundance of MapMan Categories",
              titlefont = list(size=fontsizes[3]),
              legend = list(font = list(size=fontsizes[4])))

    if(logscale==TRUE) {
      p<-layout(p,
               margin = list(b = 160),
               separators = '.,',
               xaxis=list(tickangle = 45,titlefont=list(size=fontsizes[1]),tickfont=list(size=fontsizes[2])),
              yaxis = list(title = "",type = "log",tickformat="e,",titlefont=list(size=fontsizes[1]),tickfont=list(size=fontsizes[2])),
              title = "Abundance of MapMan Categories",
              titlefont = list(size=fontsizes[3]),
              legend = list(font = list(size=fontsizes[4])))
      return(p)
    }

    return(p)
  }

  if(plottype=="pie"){
    p <- plot_ly() %>%
      add_pie(mapman,
              labels = ~rownames(mapman),
              values  = ~mapman[,1],
              name = colnames(mapman)[1],
              marker=list(colors = paleta_mapman),
              domain = list(x = c(0, 0.4), y = c(0.4, 1)),
              sort=FALSE) %>%
      add_pie(mapman,
              labels = ~rownames(mapman),
              values = ~mapman[,2],
              name = colnames(mapman)[2],
              marker=list(colors = paleta_mapman),
              domain = list(x = c(0.6, 1), y = c(0.4, 1)),
              sort=FALSE)

    return(p)
  }

  if(plottype=="heatmap"){
    require(pheatmap)
    require(RColorBrewer)
    values <- listaMMO$mean
    p<-pheatmap(values, color = colorRampPalette(rev(brewer.pal(n = 11, name = "Spectral")))(100),
                clustering_distance_rows=hmparameters[1],
                clustering_distance_cols=hmparameters[2],
                clustering_method=hmparameters[3],
                display_numbers=hmparameters[4],
                scale=hmparameters[5])
    class(p)<-"Plot_mapman_heatmap"
    return(p) #ojo, para plotearlo a antojo posteriormente usar grid::grid.draw(p$gtable)
  }
}
