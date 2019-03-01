# Aqui hay un bug. Si para un vennplot de 3 tratamientos solo hay 3 grupos (ejemplo A,B,C; BC; C), en vez de los 7 que tocan peta.
# Arreglar para ¿poner ceros al resto?. ¿Otra alternativa? Se creo un dataset bugvenn con el set de datos necesario para replicar esto.
# gruposvenn_met <- Venn_group(datosChlamyPreprocesados2,initialrow = 1,initialcolumn = 5,treatment1col = 1,treatment2col = 1,treatment = 1,omiclevel="Metabolites",threshold=2)
# vennplot_met<-Venn_plot(gruposvenn_met, alpha=0.5, num.cex=2 ,font.cex=2)
#
# Set de datos en el directorio raiz. Aprovechar y definir tests unitarios para esta funcion, empleando todas las posibilidades.
#
# Cambiar nombre de argumentos num.cex font.cex para homogeneizarlos con el resto de funciones graficas

Venn_group<-function(datalist,initialrow=1, initialcolumn=2, treatment1col=1,treatment2col=2, treatment=1,omiclevel=NULL,threshold=2){
  require(plyr)
  options(warn=-1)
  if(class(datalist)!="POL") stop("A POL object is required. Please run pre-processing steps")
  datasetnames<-names(datalist)
  
  #LAU: Esto esta muy bien pero obligas a que la matriz tenga una estructura bastante concreta, habra que tenerlo en cuenta en la documentacion, 
  # se podria hacer mas flexible pero de momento impones que los tratamientos (maximo 2 clasificaciones y una tercera columna suma de ambos) 
  #estan al principio de la tabla (/POR COjONES/) en las posiciones 1,2,3 luego lo que sea y despues ya los valores numericos
 
  #Generamos un vector para hacer el split por tratamientos
  if(treatment==1) vectortratamientos <- as.factor(datalist[[1]][initialrow:nrow(datalist[[1]]),treatment1col])
  if(treatment==2) vectortratamientos <- as.factor(datalist[[1]][initialrow:nrow(datalist[[1]]),treatment2col])
  if(treatment==3) vectortratamientos <- as.factor(paste(datalist[[1]][initialrow:nrow(datalist[[1]]),treatment1col],datalist[[1]][initialrow:nrow(datalist[[1]]),treatment2col], sep="&"))
  
  #Extraemos los valores significativos de las tablas
  datalist2 <- lapply(datalist, function(x) x[initialrow:nrow(x),initialcolumn:ncol(x)])
  if(is.null(omiclevel)==TRUE) omiclevel="all"
  #LAU:ahi te falta una llave para el if, I think
  
  # Seleccionamos los niveles omicos datos
 
  if(omiclevel=="all") {
    variablenames <- lapply(datalist2, function(x) colnames(x))
    variablenames <-unlist(variablenames)
    listaintermedia <-do.call(cbind,datalist2)  #### Esto es inconsistente. Cambia nombres de columnas, a veces, nombretabla.nombrevar, otras veces solo nombrevar
    colnames(listaintermedia) <- variablenames #### para arreglar lode arriba entre tanto.
    listaintermedia<-RemoveEmptyColumns(listaintermedia,1,1)#Esto es para evitar cosas raras.
  }  
  if(omiclevel!="all") {
    sapply(omiclevel, function(x) if(x %in% datasetnames==FALSE) stop("Error: Please select a valid datasetname"))
    
    if(length(omiclevel)>1) {  
    listaintermedia <-do.call(cbind,datalist2[omiclevel])
    variablenames <- lapply(datalist2[omiclevel], function(x) colnames(x))
    variablenames <-unlist(variablenames)
    colnames(listaintermedia) <- variablenames
    listaintermedia<-RemoveEmptyColumns(listaintermedia,1,1)#Esto es para evitar cosas raras.
    }
    if(length(omiclevel)==1){
      listaintermedia <-datalist2[[omiclevel]]
      variablenames <-colnames(listaintermedia)
      listaintermedia<-RemoveEmptyColumns(listaintermedia,1,1)#Esto es para evitar cosas raras.
    }
  }
  
  datos <- split(as.data.frame(listaintermedia),vectortratamientos,drop=F)
  
  Treat.data <- lapply(datos,function(x) apply(x,2,function(y) length(which(y != 0))))
  
  Treat.data <- Treat.data[unique(vectortratamientos)]
  presence.matrix <- (do.call("cbind", Treat.data)) 
  
  #hacemos una matriz con los nombres de los tratamientos que estan presentes, t.q. xij>=threshold-->nombredeltreatment; xij<threshold-->NA
  for(i in 1:nrow(presence.matrix)){
    for(j in 1:ncol(presence.matrix)){
      if(presence.matrix[i,j]>=threshold){
        presence.matrix[i,j]<-colnames(presence.matrix)[j]
      } 
      else {
        presence.matrix[i,j]<-NA
      }
    }
  }  
  presence.matrix<-as.data.frame(presence.matrix)
 
  #aqui hace falta algo que elimine las que no tienen abundacia suficiente para pasar el threshold en ninguno de los tratamientos
  #y un  mensajito que te avise de hay variables que son muy poco abundantes y que no se consideran como pertenecientes a nada
  if (length(which(apply(presence.matrix,1,function(x)all(is.na(x)))))>0) {
    cat("\n",length(which(apply(presence.matrix,1,function(x)all(is.na(x))))),"variables have been removed from the analysis, as were not in at least of",threshold,"replicates in any treatment.")
    presence.matrix<-presence.matrix[-which(apply(presence.matrix,1,function(x)all(is.na(x)))),]
  }else{
    cat("\nAll the variales were found to be in at least",threshold,"replicates in one or more treatments.")
  }
  grouping <- as.data.frame(apply(presence.matrix,1,function(x) auxiliarpaste_venn(x)))
  colnames(grouping)<-"groups"
  grouping<-cbind(row.names(grouping),grouping) #sip, esto es una cutrez, pero lo necesito para que dlply no me descojone la tabla y luego no haya forma de saber quien es quien
  
  Vennmatrix<-(dlply(grouping, .(groups))) #al split especial este no le sigue un unique porque me da igual el orden
 
  Vennmatrix<-do.call(rbind,Vennmatrix) #de hecho esta mejor ordenado por R que como estaba inicialmente :)

   #hasta este punto tengo una matriz de tantas observaciones como variables y cada una con su grupo asociado
  row.names(Vennmatrix)<-Vennmatrix[,1]
  finaltable<-as.data.frame(dplyr::count(Vennmatrix, Vennmatrix[,-1]))
  #finaltable nos da los grupos que hay y su frecuencia....
  #ahora lo suyo seria emparejar las variables con sus respectivos grupos...
  finallist<-sapply(finaltable[,1],function(x) row.names(grouping[grouping$groups==x,]))
  finallist <- lapply(finallist,function(x) paste(x,collapse=";"))
  Vennfinalobject<-do.call(rbind,finallist)
  Vennfinalobject<-cbind(finaltable,Vennfinalobject)
  #y asi consigo una matriz con el nombre de los grupos de Venn, que luego me podra reconocer el plot, su frecuencia y las variables que lo componen
  
  vennresult<-list(treatments=vectortratamientos,data=Vennfinalobject)
  class(vennresult) <- "vennanalysis"
  return(vennresult)
}

auxiliarpaste_venn<-function(v){  #no hay una forma directa de eliminar los NAs del string he intentado na.rm=TRUE y no va...
  aux<-c()
  for(j in 1:length(v))
    if(is.na(v[j])==FALSE){
      aux<-paste(aux,v[j],sep = "&")
    } else{
      aux<-aux
    }
  return(substring(aux, 2))
}

Venn_plot<-function(datalist, Euler.dist = FALSE, alpha=0.5, num.cex=3,font.cex=3){
  #LAU: he integrado los dos valores que me pedias, de tama?os de los n?meros y categorias, los meto en venn_plot porque alpha esta ahi tambien 
  #LAU: meto los diagramas de dos conjuntos
  
  require(VennDiagram)
  
  if(class(datalist)!="vennanalysis") stop("A Venn Analysis object is required. Please run Venn_group step first")
  treatments <-datalist[[1]]
  Vennfinalobject <-datalist[[2]]
  
      #cuantos tratamients hay?
  Treat.No<-length(unique(treatments))
  
  if(Treat.No==2){
    require(VennDiagram)
    #a<-as.vector.factor(unique(matriz[,treatments]))
    a<-as.vector(unique(matriz[,treatments]))
    #double treatment groups
    if(length(which(Vennfinalobject[,1]==paste(a[1],a[2],sep="&")))==0){p12<-0}
    else{p12<-Vennfinalobject[which(Vennfinalobject[,1]==paste(a[1],a[2],sep="&")),2]}
    
    #single treatment groups
    if(length(which(Vennfinalobject[,1]==a[1]))==0){area1<-0+p12}
    else{area1<-Vennfinalobject[which(Vennfinalobject[,1]==a[1]),2]+p12} 
    
    if(length(which(Vennfinalobject[,1]==a[2]))==0){area2<-0+p12}
    else{area2<-Vennfinalobject[which(Vennfinalobject[,1]==a[2]),2]+p12}
    
   
    if(dev.cur()>1) dev.off()
    draw.pairwise.venn(area1, area2, p12, category =
                         c(a[1],a[2]), euler.d =
                         Euler.dist, scaled = Euler.dist, inverted=FALSE,ext.text = TRUE, 
                       ext.percent = rep(0.05, 3), lwd = rep(2, 2), lty =
                         rep("blank", 2), col = rep("black", 2), fill = c("blue","green"),
                       alpha = rep(alpha, 2), label.col=rep("black",3),cex
                       = rep(num.cex,3), fontface = rep("plain", 3), fontfamily =
                         rep("serif", 3), cat.pos = c(0, 0), cat.dist =
                         c(0.05, 0.05), cat.col = c("blue","green"),
                       cat.cex = rep(font.cex, 2), cat.fontface = rep("plain", 2),
                       cat.fontfamily = rep("serif", 2), cat.just =
                         rep(list(c(0.5, 0.5)),2), cat.default.pos
                       = "outer", cat.prompts = FALSE, rotation.degree = 0,
                       rotation.centre = c(0.5, 0.5), ind = TRUE, sep.dist =
                         0.05, offset = 0, cex.prop = NULL, print.mode = "raw",
                       sigdigs = 1)
    class(myplot) <-"Venngrid"
    return(myplot)
  }
  
  else if(Treat.No==3){
    require(VennDiagram)
    a<-as.vector(unique(treatments))

    #triple treatment groups
    if(length(which(Vennfinalobject[,1]==paste(a[1],a[2],a[3],sep="&")))==0){n123<-0}
    else{n123<-Vennfinalobject[which(Vennfinalobject[,1]==paste(a[1],a[2],a[3],sep="&")),2]  }
    
    
    #double treatment groups
    if(length(which(Vennfinalobject[,1]==paste(a[1],a[2],sep="&")))==0){n12<-0+n123}
    else{n12<-Vennfinalobject[which(Vennfinalobject[,1]==paste(a[1],a[2],sep="&")),2]+n123}
    
    if(length(which(Vennfinalobject[,1]==paste(a[2],a[3],sep="&")))==0){n23<-0+n123}
    else{n23<-Vennfinalobject[which(Vennfinalobject[,1]==paste(a[2],a[3],sep="&")),2]+n123} 
    
    if(length(which(Vennfinalobject[,1]==paste(a[1],a[3],sep="&")))==0){n13<-0+n123}
    else{n13<-Vennfinalobject[which(Vennfinalobject[,1]==paste(a[1],a[3],sep="&")),2]+n123} 
    
    #single treatment groups
    if(length(which(Vennfinalobject[,1]==a[1]))==0){area1<-0+n12+n13-n123}
    else{area1<-Vennfinalobject[which(Vennfinalobject[,1]==a[1]),2]+n12+n13-n123}
    
    if(length(which(Vennfinalobject[,1]==a[2]))==0){area2<-0+n12+n23-n123}
    else{area2<-Vennfinalobject[which(Vennfinalobject[,1]==a[2]),2]+n12+n23-n123}
    
    if(length(which(Vennfinalobject[,1]==a[3]))==0){area3<-0+n13+n23-n123}
    else{area3<-Vennfinalobject[which(Vennfinalobject[,1]==a[3]),2]+n13+n23-n123} 
    
    #he aqui el plot... basico
    if(dev.cur()>1) dev.off()
    myplot<-draw.triple.venn(area1, area2, area3, n12, n23, n13, n123, category =
                       c(a[1],a[2],a[3]), rotation = 1, reverse = FALSE, euler.d =
                         Euler.dist, scaled = Euler.dist, lwd = rep(2, 3), lty =
                       rep("blank", 3), col = rep("black", 3), fill = paleta_tratamientos[1:3],
                     alpha = rep(alpha, 3), label.col = rep("black", 7), cex
                     = rep(num.cex, 7), fontface = rep("plain", 7), fontfamily =
                       rep("sans", 7), cat.pos = c(-40, 40, 180), cat.dist =
                       c(0.05, 0.05, 0.025), cat.col = paleta_tratamientos[1:3],
                     cat.cex = rep(font.cex, 3), cat.fontface = rep("plain", 3),
                     cat.fontfamily = rep("sans", 3), cat.just =
                       list(c(0.5, 1), c(0.5, 1), c(0.5, 0)), cat.default.pos
                     = "outer", cat.prompts = FALSE, rotation.degree = 0,
                     rotation.centre = c(0.5, 0.5), ind = TRUE, sep.dist =
                       0.05, offset = 0, cex.prop = NULL, print.mode = "raw",
                     sigdigs = 1)
    
    class(myplot) <-"Venngrid"
    return(myplot)
  }
  else if(Treat.No==4){
    require(VennDiagram)
    aa<-as.vector(unique(treatments))
    #quadruple treatment groups
    if(length(Vennfinalobject[which(Vennfinalobject[,1]==paste(aa[1],aa[2],aa[3],aa[4])))==0){m1234<-0}
    else{m1234<-Vennfinalobject[which(Vennfinalobject[,1]==paste(aa[1],aa[2],aa[3],aa[4],sep="&")),2]}  
    
    #triple treatment groups
    if(length(Vennfinalobject[which(Vennfinalobject[,1]==paste(aa[1],aa[2],aa[3])))==0){m123<-0+m1234}
    else{m123<-Vennfinalobject[which(Vennfinalobject[,1]==paste(aa[1],aa[2],aa[3],sep="&")),2]+m1234}
    
    if(length(Vennfinalobject[which(Vennfinalobject[,1]==paste(aa[1],aa[3],aa[4])))==0){m134<-0+m1234}
    else{m134<-Vennfinalobject[which(Vennfinalobject[,1]==paste(aa[1],aa[3],aa[4],sep="&")),2]+m1234}
    
    if(length(Vennfinalobject[which(Vennfinalobject[,1]==paste(aa[1],aa[2],aa[4])))==0){m124<-0+m1234}
    else{m124<-Vennfinalobject[which(Vennfinalobject[,1]==paste(aa[1],aa[2],aa[4],sep="&")),2]+m1234}
    
    if(length(Vennfinalobject[which(Vennfinalobject[,1]==paste(aa[2],aa[3],aa[4])))==0){m234<-0+m1234}
    else{m234<-Vennfinalobject[which(Vennfinalobject[,1]==paste(aa[2],aa[3],aa[4],sep="&")),2]+m1234}  
    
    #double treatment groups
    if(length(Vennfinalobject[which(Vennfinalobject[,1]==paste(aa[1],aa[2])))==0){m12<-0+m123+m124-m1234 }
    else{ m12<-Vennfinalobject[which(Vennfinalobject[,1]==paste(aa[1],aa[2],sep="&")),2]+m123+m124-m1234 }
    
    if(length(Vennfinalobject[which(Vennfinalobject[,1]==paste(aa[1],aa[3])))==0){m13<-0+m123+m134-m1234 }
    else{ m13<-Vennfinalobject[which(Vennfinalobject[,1]==paste(aa[1],aa[3],sep="&")),2]+m123+m134-m1234 }
    
    if(length(Vennfinalobject[which(Vennfinalobject[,1]==paste(aa[1],aa[4])))==0){m14<-0+m124+m134-m1234 }
    else{ m14<-Vennfinalobject[which(Vennfinalobject[,1]==paste(aa[1],aa[4],sep="&")),2]+m124+m134-m1234 }
    
    if(length(Vennfinalobject[which(Vennfinalobject[,1]==paste(aa[2],aa[3])))==0){m23<-0+m123+m234-m1234 }
    else{ m23<-Vennfinalobject[which(Vennfinalobject[,1]==paste(aa[2],aa[3],sep="&")),2]+m123+m234-m1234 }
    
    if(length(Vennfinalobject[which(Vennfinalobject[,1]==paste(aa[2],aa[4])))==0){m24<-0+m123+m234-m1234 }
    else{ m24<-Vennfinalobject[which(Vennfinalobject[,1]==paste(aa[2],aa[4],sep="&")),2]+m123+m234-m1234 }
    
    if(length(Vennfinalobject[which(Vennfinalobject[,1]==paste(aa[3],aa[4])))==0){m34<-0+m134+m234-m1234 }
    else{ m34<-Vennfinalobject[which(Vennfinalobject[,1]==paste(aa[3],aa[4],sep="&")),2]+m134+m234-m1234 }
    
    
    #single treatment groups
    if(length(Vennfinalobject[which(Vennfinalobject[,1]==aa[1]))==0){aream1<-0+m12+m13+m14-m123-m124-m134+m1234 }
    else{ aream1<-Vennfinalobject[which(Vennfinalobject[,1]==aa[1]),2]+m12+m13+m14-m123-m124-m134+m1234 }
    
    if(length(Vennfinalobject[which(Vennfinalobject[,1]==aa[2]))==0){aream2<-0+m12+m13+m14-m123-m124-m134+m1234 }
    else{ aream2<-Vennfinalobject[which(Vennfinalobject[,1]==aa[2]),2]+m12+m13+m14-m123-m124-m134+m1234 }
    
    if(length(Vennfinalobject[which(Vennfinalobject[,1]==aa[3]))==0){aream3<-0+m13+m23+m34-m123-m134-m234+m1234 }
    else{ aream3<-Vennfinalobject[which(Vennfinalobject[,1]==aa[3]),2]+m13+m23+m34-m123-m134-m234+m1234 }
    
    if(length(Vennfinalobject[which(Vennfinalobject[,1]==aa[4]))==0){aream4<-0+m14+m24+m34-m124-m134-m234+m1234 }
    else{ aream4<-Vennfinalobject[which(Vennfinalobject[,1]==aa[4]),2]+m14+m24+m34-m124-m134-m234+m1234 }
    
    
    if(dev.cur()>1) dev.off()
    myplot<-draw.quad.venn(aream1, aream2, aream3, aream4, m12, m13, m14, m23, m24,
                   m34, m123, m124, m134, m234, m1234, category = aa,
                   lwd = rep(2, 4), lty = rep("blank", 4), col =
                   rep("black", 4), fill = paleta_tratamientos[1:4], 
                   alpha = rep(alpha, 4),
                   label.col = rep("black", 15), cex = rep(num.cex, 15),
                   fontface = rep("plain", 15), fontfamily = rep("sans",
                   15), cat.pos = c(-15, 15, 0, 0), cat.dist = c(0.22,
                   0.22, 0.11, 0.11), cat.col = paleta_tratamientos[1:4], cat.cex
                   = rep(font.cex, 4), cat.fontface = rep("plain", 4),
                   cat.fontfamily = rep("sans", 4), cat.just =
                   rep(list(c(0.5, 0.5)), 4), rotation.degree = 0,
                   rotation.centre = c(0.5, 0.5), ind = TRUE, cex.prop =
                   NULL, print.mode = "raw", sigdigs = 3, direct.area =
                   FALSE, area.vector = 0)
    class(myplot) <-"Venngrid"
    return(myplot)
  }
  else if(Treat.No>=5){
    sections<-dim(Vennfinalobject)[1]
    #ahora genero un objeto, que sera el tama?o,es decir, la frecuencia de cada section
    sizesection<-Vennfinalobject[2]
    #lo siguiente es generar el objeto que me diga que tratamiento esta en cada sector, un boolean T/F, que luego se convirtio en una matriz de 0 y 1
    #va a tener que ser un bucle
    Treats<-as.vector.factor(unique(treatments))
    present.object<-c()
    
    for(k in 1:(length(Treats))){
      is.there<-c()
      is.present<-c()
      for(i in 1:sections){
        is.there<-grepl(Treats[k],Vennfinalobject[i,1])
        if(is.there==T){
          is.there<-1
        }else{
          is.there<-0
        }
        is.there<-is.there*sizesection[i,1]
        is.present<-rbind(is.present,is.there)
      }
      present.object<-cbind(present.object,is.present)
    }
    colnames(present.object)<-Treats
    rownames(present.object)<-Vennfinalobject[,1]
    present.object<-as.data.frame(present.object)
    
    #hasta aqui mi objeto para crear el plot
    #pues no, lo que varia son los colores de una seccion a otra, o blanco o el color que toque, asi que hago la matriz de colores
    colourpantone<-c("orange","gold","yellow","green","green4","blue","purple","red")
    #no tengo capaciad para elegir colores... asique no me voy a parar aqui, pongo los nombres mas tipicos y digo yo que R sabra cuales son
    #### HELP LUIS/MONICA!
    coloursvec<-colourpantone[1:Treat.No]
    colourmat<-c()
    for(k in 1:(length(Treats))){
      is.on<-c()
      is.off<-c()
      for(i in 1:sections){
        if(present.object[i,k]==0){
          is.on<-"white"
        }else{
          is.on<-coloursvec[k]
        }
        is.off<-rbind(is.off,is.on)
      }
      colourmat<-cbind(colourmat,is.off)
      
    }
    colourmat<-as.data.frame(colourmat)
    colnames(colourmat)<-Treats
    rownames(colourmat)<-Vennfinalobject[,1]
    #ahora si, crear tantos circulos como length(Treats)
    #install.packages('plotrix')
    require('plotrix')
    iniR<-0.2
    xsect<-as.vector(sizesection[,1])
    #finally: the multilayerplot!
    pie(0,0,x=xsect,edges=dim(matriz)[2], radius=(Treat.No+1)*iniR, col=as.vector(colourmat[,Treat.No]), border = NA,labels = '')
    for(i in ((length(Treats))-1):1){
              floating.pie(0,0,x=xsect,edges=dim(matriz)[2], radius=(i+1)*iniR, col=as.vector(colourmat[,i]), border = NA)
    }   
    floating.pie(0,0,x=1,edges=1000, radius=iniR, col="white", border = NA)
    
    legend(0.5, 5*iniR, Treats, col=as.character(coloursvec), pch=19,bty='n', ncol=2)
    
  }
    
}
