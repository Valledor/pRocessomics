
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



#' @name export_plot
#' @title A function for exporting pRocessomics generated plots to pdf format
#' @description This function exports plots to pdf format
#' @usage export_plot(plot_object, filename = "myplot.pdf", 
#' width=NULL, height=NULL, cutoff=NULL)
#' @param plot_object name of the plot in the RStudio environment to be exported
#' @param filename quoted name of the exported plot with the extension .pdf
#' @param width width of the pdf object
#' @param height height of the pdf object
#' @param cutoff numeric, cutoff for the network, only used if the plot to be exported is a network
#' @author Luis Valledor and Laura Lamelas
#' @export

#' @importFrom grDevices dev.cur dev.off pdf
#' @importFrom grid grid.draw
#' @importFrom ggplot2 ggsave
#' @importFrom mixOmics network circosPlot cimDiablo
#' @importFrom igraph write.graph
#' @importFrom plotly orca


export_plot <- function(plot_object, filename="myplot.pdf",width=NULL,height=NULL,cutoff=NULL){
  if(grDevices::dev.cur()>1) grDevices::dev.off()
  
  # if("DA" %in% class(plot_object)){
  #   if(is.null(width)==TRUE) width=8
  #   if(is.null(height)==TRUE) height=7
  #   if(filename=="myplot.pdf") filename="myplot_DAnetwork"
  #   texts_wizard("\nExporting network plot. It may take a while\n")
  #   set.seed(5881)
  #   exportnet<-mixOmics::network(plot_object, cutoff = cutoff, comp = c(1, 2))
  #   texts_wizard("Exporting gml network for Cytoscape. It may take a while\n")
  #   igraph::write.graph(exportnet$gR, file = paste(filename,".gml",sep=""), format = "gml")
  #   return()
  # }

  if(inherits(plot_object,"mixo_spls")){
    if(is.null(width)==TRUE) width=8
    if(is.null(height)==TRUE) height=7
    if(filename=="myplot.pdf") filename="myplot_DAnetwork"
    texts_wizard("\nExporting network plot. It may take a while\n")
    set.seed(5881)
    exportnet<-mixOmics::network(plot_object, cutoff = cutoff, comp = c(1, 2))
    texts_wizard("Exporting gml network for Cytoscape. It may take a while\n")
    igraph::write.graph(exportnet$gR, file = paste(filename,".gml",sep=""), format = "gml")
    return()
  }
  
  if(inherits(plot_object,"Plot_mapman_heatmap")){
    if(is.null(width)==TRUE) width=8
    if(is.null(height)==TRUE) height=7
    grDevices::pdf(file=filename,width,height,onefile=FALSE)
    grid::grid.draw(plot_object$gtable)
    grDevices::dev.off()
    texts_wizard("\nExporting Mapman-Heatmap plot. Complete\n")
    return()
  }
  
  if(inherits(plot_object,"Plot_proteins_heatmap")) {
    if(is.null(width)==TRUE) width=8
    if(is.null(height)==TRUE) height=7
    grDevices::pdf(file=filename,width,height,onefile=FALSE)
    grid::grid.draw(plot_object$gtable)
    grDevices::dev.off()
    return()
  }
  
  if(inherits(plot_object,"sunburst")){
    if(is.null(width)==TRUE) width=10
    if(is.null(height)==TRUE) height=8.5
    grDevices::pdf(file=filename,width,height,onefile=FALSE)
    Venn_plot(datalist = plot_object[[2]])
    #class(plot_object)<-"gList"
    #grid::grid.draw(plot_object)
    grDevices::dev.off()
    return()
  }
 
  if(inherits(plot_object,"Venngrid")){
    if(is.null(width)==TRUE) width=8
    if(is.null(height)==TRUE) height=7
    grDevices::pdf(file=filename,width,height,onefile=FALSE)
    class(plot_object)<-"gList"
    grid::grid.draw(plot_object)
    grDevices::dev.off()
    return()
  }
  
  if(inherits(plot_object,"singlekmeansplot")){
    if(is.null(width)==TRUE) width=8
    if(is.null(height)==TRUE) height=7
    if(filename=="myplot.pdf") filename="myplot_kmeans.pdf"
    if(width<=10) width=10
    if(height<=10) height=10
    ggplot2::ggsave(filename=filename,plot = plot_object[[2]],width=width,height=height,device="pdf")
    return()
  }
  
  if(inherits(plot_object,"kmeansmultipleplot")){
    if(is.null(width)==TRUE) width=8
    if(is.null(height)==TRUE) height=7
    if(filename=="myplot.pdf") filename="myplot_multiplekmeans.pdf"
    if(width<=10) width=10
    if(height<=10) height=10
    grDevices::pdf(file=filename,width,height,onefile=TRUE)
    texts_wizard("Exporting k-means plots. It may take a while\n")
    invisible(lapply(plot_object[[2]], print))
    grDevices::dev.off()
    return()
  }
  
  if(inherits(plot_object,"sgccda")){
    if(is.null(width)==TRUE) width=8
    if(is.null(height)==TRUE) height=7
    if(filename=="myplot.pdf") filename="myplot_network"
    texts_wizard("\nExporting network plot. It may take a while\n")
    set.seed(5881)
    exportnet<-mixOmics::network(plot_object,cex.node.name = 0.5,cutoff=cutoff,name.save=filename, blocks = c(1:(length(plot_object$names$blocks)-1)))
    texts_wizard("Exporting gml network for Cytoscape. It may take a while\n")
    igraph::write.graph(exportnet$gR, file = paste(filename,".gml",sep=""), format = "gml")
    return()
  }
  
  if(inherits(plot_object,"circosplot")){
    if(is.null(width)==TRUE) width=8
    if(is.null(height)==TRUE) height=7
    if(filename=="myplot.pdf") filename="myplot_circos.pdf"
    if(width<=10) width=10
    if(height<=10) height=10
    grDevices::pdf(file=filename,width,height)
    invisible(mixOmics::circosPlot(plot_object,cutoff=cutoff,color.blocks = plot_object$lvl.cols))
    grDevices::dev.off()
    return()
  }
  
  if(inherits(plot_object,"cimdiablo")){
    if(is.null(width)==TRUE) width=8
    if(is.null(height)==TRUE) height=7
    if(filename=="myplot.pdf") filename="myplot_cim.pdf"
    if(width<=10) width=15
    if(height<=10) height=15
    if(grDevices::dev.cur()>1) grDevices::dev.off()
    grDevices::pdf(file=filename,width,height,onefile=FALSE)
    #mixOmics::cimDiablo(plot_object,legend.position = 0,size.legend=1,margins = c(8,8))
    mixOmics::cimDiablo(
      plot_object,
      color.blocks = plot_object$lvl.cols,
      color.Y = plot_object$tr.cols,
      legend.position = 0,
      size.legend = 1,
      margins = c(8, 8)
    )
    grDevices::dev.off()
    return()
  }
  
  if(inherits(plot_object,"cim")){
    if(is.null(width)==TRUE) width=8
    if(is.null(height)==TRUE) height=7
    if(filename=="myplot.pdf") filename="myplot_cim.pdf"
    if(width<=10) width=10
    if(height<=10) height=10
    if(grDevices::dev.cur()>1) grDevices::dev.off()
    grDevices::pdf(file=filename,width,height,onefile=FALSE)
    mixOmics::cimDiablo(plot_object,legend.position = 0,size.legend=1,margins = c(8,8))
    grDevices::dev.off()
    return()
  }
  
  if(inherits(plot_object,"mapman.pie.plot")){
      if(is.null(width)==TRUE) width=1000
      if(is.null(height)==TRUE) height=750
      index.plot<- length(plot_object)
      if(is.null(filename)) filename <- "mapman_pie_plot.pdf"
      for(i in 1:index.plot){
             plotly::orca(p = plot_object[[i]], file = paste(i,"_",filename, sep=""), format="pdf",width=width, height=height,scale=3)}
  }
  
  else{
    
    if(is.null(width)==TRUE) width=1000
    if(is.null(height)==TRUE) height=750
    plotly::orca(p = plot_object, file = filename, format="pdf",width=width, height=height,scale=3)
  }
}

