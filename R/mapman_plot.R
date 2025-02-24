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

#' @name mapman_plot
#' @title Generate plots based on Mapman categories and groups (mapman_group ())
#' @description A function which is used to plot abundances of mapman_group object. Barplot, pie charts and heatmaps available
#' @usage mapman_plot(MMOlist, plottype = "bar", ploterror = TRUE, logscale = FALSE,
#' hmparameters = c("manhattan","manhattan","ward.D",TRUE,"row"), normbyelements = FALSE,
#' fontsizes = c(14,10,16,12))
#' @param MMOlist Mapman aggrupation table with sums of abundances according to mapman categories
#' @param plottype Character.Type of plot to be displayed, "bar", "heatmap" or "pie".
#' @param ploterror Boolean indicating whether errors within replicas should be shown
#' @param logscale Boolean indicating whether axis of bar plot should be log-scaled
#' @param hmparameters Character vector (length 5) determining main features of heatmap plot (if defined as plottype argument).
#' hmparameters <- c("distance measure used in row clustering","distance measure in column clustering","clustering method","displaynumbers","scale")
#' where distance measure used in row/col clustering: character, those methods supported by stats::dist ("correlation","euclidean","maximum","manhattan","camberra","binary","minkowski")
#' clustering method: character, those methods provided by stats::hclust ("ward.D","ward.D2","single","complete","average","mcquiry","median","centroid")
#' displaynumbers: boolean determining if numeric values should be indicated in heatmap cells.
#' scale: character, direction of abundance balancing ("row", "column","none"). If set to row or col the given numbers represent the percentage by bin or by sample respectively.
#' @param normbyelements Boolean indicating if MapMan classification abundances should be divided by the number of elements of each bin. Default FALSE
#' @param fontsizes Numeric vector, (length = 4) containing fontsizes of title, subtitle and axis labels
#' @return Plot
#' @author Luis Valledor and Laura Lamelas
#' @note Mapman annotation files could contain duplicates, importannotation function will ask you how to proceed.
#'



#' @importFrom plotly plot_ly layout add_trace add_pie
#' @importFrom pheatmap pheatmap
#' @export



mapman_plot <-
  function(MMOlist,
           plottype = "bar",
           ploterror = TRUE,
           logscale = FALSE,
           hmparameters = c("manhattan", "manhattan", "ward.D", TRUE, "row"),
           normbyelements = FALSE,
           fontsizes = c(14, 10, 16, 12)) {
    
    # Initial checks ----
    
      if(!inherits(MMOlist, c("MMOG", "MMO")))  
      stop("Not valid object. Please run mapman_group first")
    
    if (plottype %in% c("bar", "pie", "heatmap") == FALSE)
      stop("Please select a valid plot type")
    
    # Required objects ----
    values <- MMOlist$mean
    sd <- MMOlist$sd
    treatments <-
      factor(MMOlist$treatments, levels = unique(MMOlist$treatments))
    
    if (inherits(MMOlist, "MMO"))
       ploterror = FALSE
    
    if (normbyelements == TRUE) {
     
       if (inherits(MMOlist, "MMOG"))
        values <- MMOlist$mean / MMOlist$number_elements_category
    }
    # Custom aesthetics ----
    custom.aes <- Filter(function(x) "color.pRo" %in% class(get(x)), ls(envir= .GlobalEnv))
    if(length(custom.aes) != 0){
      for (i in 1:length(custom.aes)){
        new.aes <- get(custom.aes[i],envir=.GlobalEnv)
        
        switch(class(new.aes)[2], 
               "treatment" = paleta_tratamientos <- new.aes,
               "levels" = paleta_niveles <- new.aes, 
               "continuos" = paleta_continuo <- new.aes,
               "mapman" = paleta_mapman <- new.aes,
               "symbols" = simbolos <- new.aes, 
               "fontsizes" = fontsizes <- new.aes)
      }
      
    }
    
    
    # Bar plot ----
    if (plottype == "bar") {
      if (ploterror == T) {
        p <- plotly::plot_ly(
          as.data.frame(values),
          x = rownames(values),
          y = values[, 1],
          type = "bar",
          color = colnames(values)[1],
          colors = paleta_tratamientos[1:ncol(values)],
          showlegend = TRUE,
          error_y = list(
            value = sd[1],
            color = "#656a73",
            thickness = 1,
            width = 1,
            symmetric = TRUE
          )
        )
        if (ncol(values) > 1) {
          for (i in 2:ncol(values)) {
            p <- plotly::add_trace(
              p,
              x = rownames(values),
              y = values[, i],
              type = "bar",
              color = colnames(values)[i],
              error_y = list(
                value = sd[i],
                color = "#656a73",
                thickness = 1,
                symmetric = TRUE
              )
            )
          }
        }
      }
      
      if (ploterror == FALSE) {
        p <- plotly::plot_ly(
          as.data.frame(values),
          x = rownames(values),
          y = values[, 1],
          type = "bar",
          color = colnames(values)[1],
          colors = paleta_tratamientos[1:ncol(values)],
          showlegend = TRUE
        )
        if (ncol(values) > 1) {
          for (i in 2:ncol(values)) {
            p <- plotly::add_trace(
              p,
              x = rownames(values),
              y = values[, i],
              type = "bar",
              color = colnames(values)[i]
            )          
          }
        }
      }
      
      p <- plotly::layout(
        p,
        margin = list(b = 160),
        xaxis = list(
          tickangle = 45,
          titlefont = list(size = fontsizes[1]),
          tickfont = list(size = fontsizes[2])
        ),
        yaxis = list(
          title = "",
          type = "",
          tickformat = "e",
          titlefont = list(size = fontsizes[1]),
          tickfont = list(size = fontsizes[2])
        ),
        title = "Abundance of MapMan Categories",
        titlefont = list(size = fontsizes[3]),
        legend = list(font = list(size = fontsizes[4]))
      )
      
      if (logscale == TRUE) {
        p <- plotly::layout(
          p,
          margin = list(b = 160),
          separators = '.,',
          xaxis = list(
            tickangle = 45,
            titlefont = list(size = fontsizes[1]),
            tickfont = list(size = fontsizes[2])
          ),
          yaxis = list(
            title = "",
            type = "log",
            tickformat = "e,",
            titlefont = list(size = fontsizes[1]),
            tickfont = list(size = fontsizes[2])
          ),
          title = "Abundance of MapMan Categories",
          titlefont = list(size = fontsizes[3]),
          legend = list(font = list(size = fontsizes[4]))
        )
        return(p)
      }
      
      return(p)
    }
    # Pie plot ----
    if (plottype == "pie") {
      p <- list()
      for(i in 1:ncol(values)){
      mapman <- as.data.frame(values)
      p[[i]] <- plotly::plot_ly()
      p[[i]] <- plotly::add_pie(
        p = p[[i]],
        mapman,
        labels = ~ rownames(mapman),
        values  = ~ mapman[, i],
        name = colnames(mapman)[i],
        marker = list(colors = paleta_mapman),
        #domain = list(x = c(0, 0.4), y = c(0.4, 1)),
        title= paste("Mapman classification for Treatment ", colnames(mapman)[i], sep=""),
        sort = FALSE
      )
      
      }
      
      names(p) <- colnames(mapman)
      lapply(p, function(x) print(x))
      class(p) <- "mapman.pie.plot"
      return(p)
    }
    
    # Heatmap ----
    if (plottype == "heatmap") {
       if(hmparameters[5]=="none") {values <- MMOlist$mean}
       if(hmparameters[5]=="row") {values <- sweep(MMOlist$mean, 1, rowSums(MMOlist$mean), FUN="/")}
       if(hmparameters[5]=="col") {values <- sweep(MMOlist$mean, 2, colSums(MMOlist$mean), FUN="/")}
      aux_names <- MMOlist$number_bin_elements
      rownames(values) <- paste(rownames(values), "(",aux_names, ")", sep="")
      p <-
        pheatmap::pheatmap(
          values,
          color = paleta_continuo,
          clustering_distance_rows = hmparameters[1],
          clustering_distance_cols = hmparameters[2],
          clustering_method = hmparameters[3],
          display_numbers = hmparameters[4],
          angle_col = "0" #or 90 for larger names
          #scale = hmparameters[5]
        )
      class(p) <- "Plot_mapman_heatmap"
      return(p) 
    }
  }
