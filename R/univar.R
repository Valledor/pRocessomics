#
# Como estan todas las funciones fuera, no se hasta que punto merece la pena perder el tiempo en paralelizar
#
univariate <-function(datalist, initialrow=1, initialcolumn=3, treatment1col=1, 
                      treatment2col=2, treatment=1, omiclevel=NULL,parametric=TRUE, 
                      posthoc=TRUE, FDR=TRUE,round=3){
  options(warn=-1)
  #supprestartmessages o algo asi
  require(dplyr)
  require(plyr)
  require(FSA)
  
  # Checks iniciales
  if(class(datalist)!="POL") stop("\nA POL object is required. Please run pre-processing steps")
  datasetnames<-names(datalist)

  #Nos aseguramos de quitar las columnas vacias
  datalist <- lapply(datalist, function(x) RemoveEmptyColumns(x, initialrow=initialrow,initialcolumn=initialcolumn))

  #Extraemos los valores significativos de las tablas
  listadatos2 <- lapply(datalist, function(x) x[initialrow:nrow(x),initialcolumn:ncol(x)])

  #Definimos el vector que utilizaremos para comparar
  ##Lau says: entiendo que datalist[[1]][,treatment1col]==datalist[[i]][,treatment1col] and so on...
  if(treatment==1) splitvector <- as.vector(datalist[[1]][,treatment1col])
  if(treatment==2) splitvector <- as.vector(datalist[[1]][,treatment2col])
  if(treatment==3) splitvector <- as.vector(paste(datalist[[1]][,treatment1col], datalist[[1]][,treatment2col],sep="&"))
  #print(splitvector)
  #Meter aqui: if(parametric==NULL){parametric<-lapply(datalist, function(x) areUparametric())}
  
  # A partir de aqui para cada lista definiremos lo siguiente
  stattest <- c()
  if(length(unique(splitvector))==1) stop("Two or more treatments are required")
  if(length(unique(splitvector))==2) stattest<-"t_or_u"
  if(length(unique(splitvector))>=3) stattest<-"anova_or_kruskal"
  #print(stattest)
  if(stattest=="t_or_u"){
    if(parametric==TRUE){
      pvaluelist<-lapply(datalist, function(x) testT(x,splitvector,initialcolumn,initialrow))  
      if(FDR==TRUE) pvaluelist <- lapply(results, function(x) qvalues(x, round))
    }
    if(parametric==FALSE){
      pvaluelist<-lapply(datalist, function(x) testU(x,splitvector,initialcolumn,initialrow))  
      if(FDR==TRUE) pvaluelist <- lapply(results, function(x) qvalues(x, round))
    }
  
  #return(pvaluelist)
  }
  
  if(stattest=="anova_or_kruskal"){
    if(parametric==TRUE){
      pvaluelist<- lapply(listadatos2,function(y) apply(y,2,function(x) funcionanova(x,splitvector,posthoc,round)))
      if(FDR==TRUE) {pvaluelist <- lapply(pvaluelist, function(x) qvalues(x, round))
      }
    }
    if(parametric==FALSE) {
        pvaluelist<- lapply(listadatos2,function(y) apply(y,2,function(x) kruskalwallis(x,splitvector,posthoc,round)))
        if(FDR==TRUE) {pvaluelist <- lapply(pvaluelist, function(x) qvalues(x, round))
        }
    #return(pvaluelist)
    }
  }
  
  #tablas de medias y SDs
  meansdlist <- lapply(listadatos2,function(x) meansd(x,splitvector,round))

  #componemos tabla final
  compositelist <- meansdlist
  for(i in 1:length(compositelist)){
    if (FDR==F&posthoc==F) {
      pvaluelist[[i]]<-as.data.frame(pvaluelist[[i]])
      colnames(pvaluelist[[i]])<-"p.value"
      compositelist[[i]] <- cbind(compositelist[[i]], pvaluelist[[i]])}
    else {compositelist[[i]] <- cbind(compositelist[[i]], t(pvaluelist[[i]]))}
  }

  results <- list(meansd=meansdlist, pvalue=pvaluelist, composite=compositelist)
  class(results) <- "UNIVAR"
 return(results)
  

  return(compositelist)
}


funcionanova <- function(vector,tratamiento,posthoc=TRUE,round=5){
  res.anova=aov(vector~as.factor(tratamiento), na.rm=TRUE)
  pvalor <- anova(res.anova)$'Pr(>F)'[1:1]
  names(pvalor)[1] <- "p-value"
  valordevuelta <-as.matrix(pvalor)
  if(posthoc==TRUE){
    res.tukey <- TukeyHSD(res.anova)
    pvalores.tukey <-res.tukey[[1]][,4]
    valordevuelta<-c(pvalor,pvalores.tukey)
    valordevuelta <- round(valordevuelta,round)
    return(valordevuelta)
  }
  valordevuelta <- t(round(valordevuelta,round))
  return(valordevuelta)
}

kruskalwallis <- function(vector,tratamiento,posthoc=TRUE,round=5){
  reskruskal<-kruskal.test(vector~as.factor(tratamiento))
  pvaloresKruskal<-reskruskal$p.value
  names(pvaloresKruskal)[1] <- "p-value"
  valordevuelta <-as.matrix(pvaloresKruskal)
  if(posthoc==TRUE){
    require("FSA")
    dt<-dunnTest(vector~as.factor(tratamiento),method = "bonferroni")
    pvaldunn<-dt[[2]][,3]
    names(pvaldunn)<-dt[[2]][,1]
    valordevuelta<-c(pvaloresKruskal,pvaldunn)
    valordevuelta <- round(valordevuelta,round)
    return(valordevuelta)
  }
  
  valordevuelta <- t(round(valordevuelta,round))
  return(valordevuelta)
}

qvalues <- function(x,round=5){
  # This transform vector in data frame, when only p values and not Kruskal are requested
  if(class(x)=="numeric"){
    x <-t(as.data.frame(x))
    row.names(x) <- "p-value"
  }
  qval<-p.adjust(x[1,],method = "BH")
  qval<-round(qval,round)
  if(nrow(x)==1) y <- rbind(x[1,],qval)
  if(nrow(x)>=2) y <- rbind(x[1,],qval,x[2:nrow(x),])
  row.names(y)[1:2] <- c(rownames(x)[1],"q-value")
  return(y)
}

meansd <- function(table,splitvector,decimals=3){
  listpertreatments <- split(table,splitvector,drop=T)
  means<-lapply(listpertreatments, colMeans)
  means<-as.data.frame(do.call(cbind, means))
  sds <- lapply(listpertreatments, function(x) apply(x,2,sd))
  sds<-as.data.frame(do.call(cbind,sds))
  means<-round(means,decimals)
  sds<-round(sds,decimals)
  plusminus <- rep(c("±"),nrow(means))
  outputtable<-c()
  for(i in 1:ncol(means)){
    temptable<-cbind(means[,i],plusminus,sds[,i])
    colnames(temptable)<- c(paste("Mean-",(colnames(means)[i]),sep="")," ",paste("SD-",(colnames(means)[i]),sep=""))
    outputtable <-cbind(outputtable,temptable)
  }
  row.names(outputtable)<-colnames(table)
  return(outputtable)
}

testT<-function(matriz,splitvector, initialcolumn, initialrow){
  #replicas por tratamiento?
  nrep<-count(splitvector)[1,2]
  #haria falta un check de que todos los tratamientos de todas las variables tienen el mismo numero de replicas
  #Loop
  p.fin.t<-c()
  for(j in initialcolumn:ncol(matriz)){
    a<-t.test(matriz[initialrow:(initialrow+nrep-1),j],matriz[initialrow+nrep:nrow(matriz),j])$p.value
    p.fin.t<-c(p.fin.t,a)
    j=j+1
  }
  p.fin.t<-as.data.frame(t(p.fin.t))
  colnames(p.fin.t)<-colnames(matriz[initialcolumn:ncol(matriz)])
  rownames(p.fin.t)<-"t.test"
  return(p.fin.t)
}

testU<-function(matriz,splitvector, initialcolumn, initialrow){
  #replicas por tratamiento?
  nrep<-count(splitvector)[1,2]
  #haria falta un check de que todos los tratamientos de todas las variables tienen el mismo numero de replicas
  #Loop
  p.fin.t<-c()
  for(j in initialcolumn:ncol(matriz)){
    a<-wilcox.test(matriz[initialrow:(initialrow+nrep-1),j],matriz[initialrow+nrep:nrow(matriz),j])$p.value
    p.fin.t<-c(p.fin.t,a)
    j=j+1
  }
  p.fin.t<-as.data.frame(t(p.fin.t))
  colnames(p.fin.t)<-colnames(matriz[initialcolumn:ncol(matriz)])
  rownames(p.fin.t)<-"Wilcox.test"
  return(p.fin.t)
}
