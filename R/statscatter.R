#' @importFrom RColorBrewer brewer.pal 
#' @importFrom plotly plot_ly colorbar layout


statscatter<-function(matriz,name,savetofile,filename,histogram=T){
  s<-as.data.frame(matriz)
  if(histogram == FALSE){
  ax <- list(
    title = "Variables",
    showgrid=F,
    showticklabels=F,
    zeroline=TRUE,
    showline=TRUE
  )
  
  ay <- list(
    title="p-values",
    range="autorange", 
    showgrid=T,
    showticklabels=T,
    zeroline="black",
    zeroline = TRUE,
    showline = TRUE
  )
  
  colorscale<-RColorBrewer::brewer.pal(11,"RdYlGn")
  if(ncol(s)==1){b<-paste(rownames(s),s$"p-value")}
  if(ncol(s)==2){b<-paste(rownames(s), '<br>',s$"DESCRIPTION", '<br>',s$"p-value")}
  p <- plotly::plot_ly(s, y = s[,1], x = rownames(s), mode='markers', type='scatter', 
                       color=s[,1], colors=colorscale,
                       hoverinfo="text",
                       hovertext= b) 
  
  p<-plotly::colorbar(p, limits = c(0.0, 1.0),which = 1) 
  
  p<-plotly::layout(p, title= paste("p value data distribution of",name),
                    shapes=list(type='line', x0= min(rownames(s)), x1= max(rownames(s)), y0=0.05, y1=0.05, line=list(color = 'black',dash='dash', width=3)),
                    xaxis = ax,
                    yaxis = ay, 
                    showlegend=F)
  
  class(p) <- c("plotly","htmlwidget","univarplot")
  
  if(savetofile==F){
    print(p)
    texts_wizard("\nPress [enter] to continue/display next plot")
    invisible(readline())
  }
  else{
    filename <-paste(sub('\\.pdf$', '', filename) , "_",name, ".pdf",sep = "")
    export_plot(plot_object = p,filename = filename )
    texts_wizard(as.character(paste("\nExporting ",filename," to ",getwd(), sep="")))
  }
  
  return(p)}
  else{
    
    p <- plotly::plot_ly(
      x = s[,1], 
      mode="markers",
      type = "histogram", 
      xbins=list(
        start = 0, end=1, size=0.05),
      marker = list(
        color=c("#a50026","#d73027","#e31a1c","#f46d43","#ff7f00",
                "#fdae61","#fee08b","#ffed6f","#ffffbf","#ffff99",
                "#e6f598","#d9ef8b","#a6d96a","#66bd63","#1a9850",
                "#276419","#006837","#00441b","#762a83","#40004b")))
    
    p <- plotly::layout(p, title= paste("p value data distribution of",name),
                        xaxis= list(title="p-value ranges"),yaxis=list(title="frequency"),showlegend=F)            
    class(p) <- c("plotly","htmlwidget","univarplot")
    if(savetofile==F){
      print(p)
      texts_wizard("\nPress [enter] to continue/display next plot")
      invisible(readline())
    }
    else{
      filename <-paste(sub('\\.pdf$', '', filename) , "_",name, ".pdf",sep = "")
      export_plot(plot_object = p,filename = filename )
      texts_wizard(as.character(paste("\nExporting ",filename," to ",getwd(), sep="")))
    }
    return(p)
   
  }
}
