export_table <- function(result_list,filename=NULL){
  library(xlsx)
  if(class(result_list)=="MMOG"){
    cat("EXPORTING MAPMAN TABLES\nA single core will be used. It may take a while...\n ")
    cat("Generating and saving Mapman table. Data was grouped. Mean and SD will be provided.\n")
    colnames(result_list$mean) <-paste(colnames(result_list$mean),"-Mean",sep="")
    colnames(result_list$sd) <-paste(colnames(result_list$sd),"-SD",sep="")
    newtable <-cbind(result_list$mean[1],result_list$sd[1])
    if(ncol(result_list$mean)>1){
      for(i in 2:ncol(result_list$mean)){
        newtable <-cbind(newtable,result_list$mean[i],result_list$sd[i])
      }
    }
    if(is.null(filename)) filename= "mapman_table.xlsx"
    write.xlsx(newtable, file=filename, sheetName="MapmanGroups")
    cat("Done!\n")
  }
  if(class(result_list)=="MMO"){
    cat("EXPORTING MAPMAN TABLES\nA single core will be used. It may take a while...\n ")
    cat("Generating and saving Mapman table. Data was ungrouped. Mean will be provided.\n")
    colnames(result_list$mean) <-paste(colnames(result_list$mean),"-Mean",sep="")
    if(is.null(filename)) filename= "mapman_table.xlsx"
    write.xlsx(result_list$mean, file=filename, sheetName="Mapman")
    cat("Done!\n")
  }
  if(class(result_list)=="pcaanalysis"){
  cat("EXPORTING PCA TABLES\nA single core will be used. It may take a while...\n ")
  cat("Generating and saving Scores table. Step 1/3\n")
  if(is.null(filename)) filename= "pca_table.xlsx"
  write.xlsx(result_list$pca$x, file=filename, sheetName="Scores")
  cat("Generating and saving Loadings table. Step 2/3\n")
  write.xlsx(result_list$pca$rotation, file=filename, sheetName="Loadings", append=TRUE)
  cat("Generating and saving Explained Variance table. Step 3/3\n")
  exportpcatable<-rbind(round((result_list$pca$sdev^2)*100/sum(result_list$pca$sdev^2),2),cumsum(round((result_list$pca$sdev^2)*100/sum(result_list$pca$sdev^2),2)),result_list$pca$sdev)
  row.names(exportpcatable)<-c("Proportion_of_Variance","Cumulative_Proportion_of_Variance","Standard_Deviation")
  colnames(exportpcatable)<-paste("PC",c(1:length(result_list$pca$sdev)),sep="")
  write.xlsx(exportpcatable, file=filename, sheetName="Proportion of Variance", append=TRUE)
  }
  if(class(result_list)=="icaanalysis"){
  cat("EXPORTING ICA TABLES\nA single core will be used. It may take a while...\n ")
  cat("Generating and saving signal estimates table. Step 1/3\n")
  if(is.null(filename)) filename= "ica_table.xlsx"
  write.xlsx(result_list$ica$S, file=filename, sheetName="Signal estimates")
  cat("Generating and saving estimated unmixing matrix. Step 2/3\n")
  write.xlsx(result_list$ica$W, file=filename, sheetName="Unmixing matrix", append=TRUE)
  cat("Generating and saving variance-accounted-for by each component. Step 3/3\n")
  exportpcatable<-data.frame(ICs=paste("IC",c(1:length(result_list$ica$vafs)),sep=""), ExplVar=round(result_list$ica$vafs*100,2), Cumulat=round(cumsum(result_list$ica$vafs*100),2))
  write.xlsx(exportpcatable, file=filename, sheetName="Variance accounted", append=TRUE)
  }
  if(class(result_list)=="mcoaanalysis"){
  cat("EXPORTING MCOA TABLES\nA single core will be used. It may take a while...\n ")
  if(is.null(filename)) filename= "mcoa_table.xlsx"
  cat("Generating and saving synthetic scores table. Step 1/5\n")
  write.xlsx(result_list$mcoa$mcoa$SynVar, file=filename, sheetName="SynVar Scores")
  cat("Generating and saving co-inertia normed scores. Step 2/5\n")
  write.xlsx(result_list$mcoa$mcoa$Tl1, file=filename, sheetName="Co-inertia Norm Score", append=TRUE)
  cat("Generating and saving column coordinates onto synthetic scores. Step 3/5\n")
  write.xlsx(result_list$mcoa$mcoa$Tco, file=filename, sheetName="Coordinates onto Syn Scores", append=TRUE)
  cat("Generating and saving pseudo-eigenvalues. Step 4/5\n")
  write.xlsx(result_list$mcoa$mcoa$cov2, file=filename, sheetName="Pseudo-eigenvalues", append=TRUE)
  norm_eigenvalues <- result_list$mcoa$mcoa$pseudoeig/sum(result_list$mcoa$mcoa$pseudoeig) #Sacamos los autovalores, y os dividimos por la suma para tener el porcentaje de varianza en teoria recogida por cada tratamiento
  screeplotdata <- data.frame(Cs=paste("Comp",c(1:length(norm_eigenvalues)),sep=" "), ExplVar=norm_eigenvalues, Cumulat=cumsum(norm_eigenvalues)*100/sum(norm_eigenvalues),stringsAsFactors = FALSE) #Creamos tabla de datos con % varianza y % var acumulado
  screeplotdata$Cs <- factor(screeplotdata$Cs, levels = unique(screeplotdata$Cs)[order(screeplotdata$ExplVar, decreasing = TRUE)]) # Ordenamos los PCs segun valor
  cat("Generating and saving cumulative pseudo-eigenvalues. Step 5/5\n")
  write.xlsx(screeplotdata , file=filename, sheetName="Cumulative pseudo-eigenvalues", append=TRUE)
  }
  if(class(result_list)=="vennanalysis"){
    cat("EXPORTING VENN GROUPING VARIABLE TABLE\nA single core will be used. It may take a while...\n ")
    if(is.null(filename)) filename= "venn_table.xlsx"
    cat("Saving xlsx file. Step 1/1\n")
    write.xlsx(result_list$data, file=filename, sheetName="Venn Groups")
  }
  if(class(result_list)=="kmeansanalysis"){
    cat("EXPORTING K-MEANS GROUPING TABLE\nA single core will be used. It may take a while...\n ")
    if(is.null(filename)) filename= "kmeans_table.xlsx"
    cat("Saving xlsx file. Step 1/1\n")
    write.xlsx(result_list$kmeans_group, file=filename, sheetName="Kmeansgroups")
}
  if(class(result_list)=="splsanalysis"){
    spls_object<-result_list[[1]]
    cat("EXPORTING SPLS TABLES\nA single core will be used. It may take a while...\n ")
    cat("Generating and saving scores tables. Step 1/3\n")
    if(is.null(filename)) filename= "spls_table.xlsx"
    write.xlsx(spls_object$variates$X, file=filename, sheetName="Scores-X")
    write.xlsx(spls_object$variates$Y, file=filename, sheetName="Scores-Y")
    cat("Generating and saving variable loadings and scores tables. Step 2/3\n")
    write.xlsx(spls_object$loadings$X, file=filename, sheetName="Loadings-X", append=TRUE)
    write.xlsx(spls_object$loadings$Y, file=filename, sheetName="Loadings-Y", append=TRUE)
    vartable <-invisible(plotVar(spls_object))
    write.xlsx(vartable[,1:3], file=filename, sheetName="VarScores", append=TRUE)
    cat("Generating and saving Explained Variance table. Step 3/3\n")
    write.xlsx(unlist(spls_object$explained_variance), file=filename, sheetName="Expl.Var", append=TRUE)
  }
  if(class(result_list)=="daanalysis"){
    spls_object<-result_list[[1]]
    cat("EXPORTING SPLS TABLES\nA single core will be used. It may take a while...\n ")
    cat("Generating and saving scores tables. Step 1/3\n")
    if(is.null(filename)) filename= "spls_table.xlsx"
    write.xlsx(spls_object$variates$X, file=filename, sheetName="Scores-X")
    write.xlsx(spls_object$variates$Y, file=filename, sheetName="Scores-Y")
    cat("Generating and saving variable loadings and scores tables. Step 2/3\n")
    write.xlsx(spls_object$loadings$X, file=filename, sheetName="Loadings-X", append=TRUE)
    write.xlsx(spls_object$loadings$Y, file=filename, sheetName="Loadings-Y", append=TRUE)
    vartable <-invisible(plotVar(spls_object))
    write.xlsx(vartable[,1:3], file=filename, sheetName="VarScores", append=TRUE)
    cat("Generating and saving Explained Variance table. Step 3/3\n")
    write.xlsx(unlist(spls_object$explained_variance), file=filename, sheetName="Expl.Var", append=TRUE)
  }
  if(class(result_list)=="diabloanalysis"){
    diablo_object<-result_list[[1]]
    cat("EXPORTING SPLS TABLES\nA single core will be used. It may take a while...\n ")
    cat("Generating and saving scores tables. Step 1/3\n")
    if(is.null(filename)) filename= "diablo_table.xlsx"
    for(i in 1:length(diablo_object$variates)) {
      write.xlsx(diablo_object$variates[[i]], file=filename, sheetName=paste("Score-",names(diablo_object$explained_variance)[i],sep=""))
    }
    cat("Generating and saving variable loadings and scores tables. Step 2/3\n")
    for(i in 1:length(diablo_object$loadings)) {
      write.xlsx(diablo_object$loadings[[i]], file=filename, sheetName=paste("Loading-",names(diablo_object$explained_variance)[i],sep=""), append=TRUE)
    }
    vartable <-invisible(plotVar(diablo_object))
    write.xlsx(vartable[,1:3], file=filename, sheetName="VarScores", append=TRUE)
    cat("Generating and saving Explained Variance table. Step 3/3\n")
    write.xlsx(unlist(diablo_object$explained_variance), file=filename, sheetName="Expl.Var", append=TRUE)
}
}

export_plot <- function(plot_object, filename="myplot.pdf",width=NULL,height=NULL,cutoff=NULL){
  if(dev.cur()>1) dev.off()
  if(is.null(width)==TRUE) width=8
  if(is.null(height)==TRUE) height=7
  if(class(plot_object)=="Plot_mapman_heatmap"){
    pdf(file=filename,width,height,onefile=FALSE)
    grid::grid.draw(plot_object$gtable)
    dev.off()
    return()
  }
  else if(class(plot_object)=="Plot_proteins_heatmap"){
      pdf(file=filename,width,height,onefile=FALSE)
      grid::grid.draw(plot_object$gtable)
      dev.off()
      return()
  }
  else if(class(plot_object)=="vennanalysis"){
    pdf(file=filename,width,height,onefile=FALSE)
    class(plot_object)<-"gList"
    grid.draw(plot_object)
    dev.off()
    return()
  }
  else if(class(plot_object)=="singlekmeansplot"){
    if(filename=="myplot.pdf") filename="myplot_kmeans.pdf"
    if(width<=10) width=10
    if(height<=10) height=10
    ggsave(filename=filename,plot = plot_object[[2]],width=width,height=height,device="pdf")
    return()
  }
  else if(class(plot_object)=="kmeansmultipleplot"){
    if(filename=="myplot.pdf") filename="myplot_multiplekmeans.pdf"
    if(width<=10) width=10
    if(height<=10) height=10
    pdf(file=filename,width,height,onefile=TRUE)
    cat("Exporting k-means plots. It may take a while\n")
    invisible(lapply(plot_object[[2]], print))
    dev.off()
    return()
  }
  else if("networkplot"%in%class(plot_object)==TRUE){
    library(igraph)
    if(filename=="myplot.pdf") filename="myplot_network"
    cat("Exporting network plot. It may take a while\n")
    set.seed(5881)
    exportnet<-network(plot_object,cex.node.name = 0.5,cutoff=cutoff,name.save=filename)
    cat("Exporting gml network for Cytoscape. It may take a while\n")
    write.graph(exportnet$gR, file = paste(filename,".gml",sep=""), format = "gml")
  }
  else if("networkplotannot"%in%class(plot_object)==TRUE){
    library(igraph)
    if(filename=="myplot.pdf") filename="myplot_annotated_network"
    cat("Exporting network plot. It may take a while\n")
    set.seed(5881)
    exportnet<-network(plot_object,cex.node.name = 0.5,cutoff=cutoff,save="pdf",name.save=filename)
    cat("Exporting gml network for Cytoscape. It may take a while\n")
    write.graph(exportnet$gR, file = paste(filename,".gml",sep=""), format = "gml")
  }
  else if("circosplot"%in%class(plot_object)){
      if(filename=="myplot.pdf") filename="myplot_circos.pdf"
      if(width<=10) width=10
      if(height<=10) height=10
      pdf(file=filename,width,height)
      invisible(circosPlot(plot_object,cutoff=cutoff))
      dev.off()
    }

  else if("cimdiablo"%in%class(plot_object)){
    if(filename=="myplot.pdf") filename="myplot_cim.pdf"
    if(width<=10) width=15
    if(height<=10) height=15
    if(dev.cur()>1) dev.off()

    pdf(file=filename,width,height,onefile=FALSE)
      cimDiablo(plot_object,legend.position = 0,size.legend=1,margins = c(8,8))
    dev.off()
  }
    else if("cim"%in%class(plot_object)){
      if(filename=="myplot.pdf") filename="myplot_cim.pdf"
      if(width<=10) width=10
      if(height<=10) height=10
      if(dev.cur()>1) dev.off()
      pdf(file=filename,width,height,onefile=FALSE)
      cimDiablo(plot_object,legend.position = 0,size.legend=1,margins = c(8,8))
      dev.off()
    return()
  }

  else{

  if (!require("webshot")) install.packages("webshot")
  library(webshot)
  if(is.null(width)==TRUE) width=992
  if(is.null(height)==TRUE) height=744
  export(plot_object, file = filename, vwidth=width, vheight=height)
  }
}

