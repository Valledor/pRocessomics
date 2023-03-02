#
# Copyright 2019 Luis Valledor / Laura Lamelas
# Last revision 01.2023
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

#' @importFrom stats ks.test pnorm
#' @importFrom car leveneTest

parametricity<-function(datalist, initialcolumn, initialrow, treatment1col, treatment2col, treatment,name=NULL){
  #suppressPackageStartupMessages(require(dplyr))
  #suppressPackageStartupMessages(require(plyr))
  #suppressPackageStartupMessages(require(car))
  options(warn = -1)
  if(!inherits(datalist, "POL")) stop("A POL object is required. Please run pre-processing steps")
  if(is.null(name)){
    datalist.name<-deparse(substitute(datalist))}
  if(!is.null(name)) {datalist.name<-name}
  datasetnames <- names(datalist)
  datalist <- lapply(datalist, function(x) RemoveEmptyColumns_univar(x, initialrow=initialrow,initialcolumn=initialcolumn))
  listadatos2 <- lapply(datalist, function(x) x[initialrow:nrow(x),initialcolumn:ncol(x)])
  
  #grouping by treatment
  if(treatment==1) splitvector <- as.vector(datalist[[1]][,treatment1col])
  if(treatment==2) splitvector <- as.vector(datalist[[1]][,treatment2col])
  if(treatment==3) splitvector <- as.vector(paste(datalist[[1]][,treatment1col], datalist[[1]][,treatment2col],sep="&"))
  splitvector<-factor(splitvector,levels=unique(splitvector))
  
  listadatos3<-lapply(listadatos2, function(x) split(x,f=splitvector))
  ks.list<-list()
  
  for(j in 1:length(listadatos3)){
    aux<-c()
    ks.res<-lapply(listadatos3[[j]], function(x)  apply(as.data.frame(x),2, function(y) as.data.frame(stats::ks.test(y,stats::pnorm)$p.value)))
    ks.res2<-lapply(ks.res, function(x) do.call(rbind,x))
    for(k in 1:length(unique(splitvector))){
      aux<-cbind(aux,ks.res2[[k]][,1])
    }
    rownames(aux)<-rownames(ks.res2[[j]])
    colnames(aux)<-paste(unique(splitvector),"KS-p.value",sep="_")
    aux<-as.data.frame(aux)
    
    ks.list[[j]]<-aux
    
    NOT.normals<-length(which(aux<=0.1))
    all<-(ncol(as.data.frame(listadatos2[[j]])))*length(unique(splitvector))
    normals<-all-NOT.normals
    texts_wizard(paste("\nTesting normality of", names(datalist)[[j]],"dataset.",sep=" "))
    texts_wizard(paste("Out of ",all," groups (treatment * variables), ", normals, " (",round(normals/all*100,2), "%) are normally distributed (Kolmogorov-Smirnov p>0.1).",sep=""))
  }
  names(ks.list)<-datasetnames
  
  levenes2<-lapply(listadatos2, function(x) apply (as.data.frame(x),2, function(y) car::leveneTest(y,splitvector)$'Pr(>F)')) 
  levenes<-lapply(levenes2, function(x) (as.data.frame(x[1,])))
  for(i in 1:length(datalist)){
    colnames(levenes[[i]])<-c("p.value_Levene")
    parametrics<-length(which(levenes[[i]][,1]>=0.05))
    all<-ncol(as.data.frame(listadatos2[[i]]))
    texts_wizard(paste("\nTesting variances equality in", names(datalist)[[i]],"dataset.",sep=" "))
    texts_wizard(paste("Out of ",all," variables, ",parametrics,  " (",round(parametrics/all*100,2), "%) have equal variances (Levene p > 0.05)",sep=""))
  }
  # homocedascity_plot_ly(levenes = levenes, datalist.name)
  plot<-parametricity_plot_ly(ks.list = ks.list,levenes = levenes,datalist.name)
  plot
  texts_wizard("\nPlot legend:\nParametricity plot depicts whether your data is (or not) normally distributed according to the different treatments and whether the variance of the data is homocedastic (light green) or heterocedastic (yellow).")
  texts_wizard("Normality is illustrated taking into account each treatment of each variable. The color cases show how many treatments of each variable are normally distributed.")
  texts_wizard("For considering a dataset is parametric each line (normality and homocedasticity) should be the greener the better!") #esto lo vas a quitar pero no se me ocurre nada mejor
  texts_wizard("\nDo you want to save this plot and associated tables?")
  saveplot<-select.list.edit(c("Yes","No"), preselect = NULL, multiple = F, title = NULL, graphics = FALSE)
  if(saveplot=="Yes"){
    texts_wizard("Saving...")
    diahora <-gsub(":","",gsub(" ","_", gsub("-","",Sys.time())))
    diahora <-substr(diahora,3,nchar(diahora)-2)
    plotfilename <- paste(diahora,"_parametricity_plot.pdf",sep="")
    export_plot(plot_object = plot, filename = plotfilename)
    res.parametricity<-list(ks.list,levenes)
    names(res.parametricity)<-c("normal.distribution", "homocedasticity")
    class(res.parametricity)<-"parametricity"
    tablefilename<-paste(diahora,"_parametricity_table.xlsx",sep="")
    export_table(result_list = res.parametricity, filename = tablefilename)
    texts_wizard(paste("\nPlot and table saved in your working directory,",paste(getwd(),", as",sep=""),plotfilename,"and", tablefilename,".", sep=" "))
    }
  
  return()
}

