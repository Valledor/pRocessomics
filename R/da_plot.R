

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

#' @name da_plot
#' @title A function to plot a Discriminant Analysis (DA)
#' @description This function allows to fully depict the results from the Discriminant analysis
#' @usage da_plot(da_list, treatment = 1, plottype = "Distance", variablespace = "x",
#' useannot = FALSE, fontsizes = c(14,10,16,12), cutoff = 0.8, fortopscoring = c(1,20,"abs"), 
#' confidence = 0.9)
#' @param da_list da analysis object, class "daanalysis"
#' @param treatment numeric value to determine the colors and symbols to be displayed in the plot, used in Score plot and Topscring plot. Available optiones are:
#' \itemize{
#' \item 1 to show different colors according to treatment 1 (treatment1col)
#' \item 2 to show different colors according to treatment 2 (treatment2col)
#' \item 3 to show different colors according to the combination of both treatments
#' \item 4 to show different colors according to treatment 1 and different symbols for treatment 2
#' }
#' @param plottype available options:
#' \itemize{
#' \item Distance, sample plot each dot represent a sample
#' \item Distance_Ellipse, distance plot with confidence ellipses by treatment
#' \item Topscoring, most representantive loadings according to Distance plot
#' \item Var, variable plot, each dot represent a variable. select dataset and annotation
#' \item Composite, combination of the two most representative plots
#' \item Composite_Annotated, same as composite but with annotations, recommended for exploratory analysis
#' \item Network, network
#' \item All all aforementioned plots will be exported as .pdf files in your working directory.
#' }
#' @param variablespace c("x","y","xy")
#' @param useannot logical, indicating if annotations should be used or not, please note the previous plsda analysis should include annotations
#' @param fontsizes font sizes for plot elements
#' @param cutoff cutoff value for network prunning
#' @param fortopscoring vector indicating the number of the component to be analyzed, the top n loadings to plot and the mode to compute "abs" for absolute value and "posneg" for positive and negative loading analyzed separately
#' @param confidence Numeric value between 0 and 1 for confidence ellipses in Distance_Ellipse plot
#' @return plotly type plot
#' @author Luis Valledor and Laura Lamelas
#' @export
#' @importFrom plotly plot_ly subplot layout style ggplotly
#' @importFrom dplyr left_join
#' @importFrom mixOmics network plotVar
#' @importFrom utils capture.output
#' @importFrom ggplot2 ggplot theme scale_color_manual scale_fill_manual stat_ellipse geom_point theme_minimal aes
#' @importFrom stats na.omit
da_plot <-
  function(da_list,
           treatment = 1,
           plottype = "Distance",
           variablespace = "x",
           useannot = FALSE,
           fontsizes = c(14, 10, 16, 12),
           cutoff = 0.8,
           fortopscoring = c(1, 20, "abs"),
           confidence = 0.9) {
    
  # Checks ----
 
  if(!inherits(da_list, "daanalysis"))  
    stop("DA analysis object is expected, please run da_analysis first")
  if (plottype %in% c("Composite", "Distance","Distance_Ellipse", "Var", "Topscoring", "Network", "All", "Composite_Annotated") == FALSE)
    stop("Please select a valid plot type (\"Composite\", \"Distance\",\"Distance_Ellipse\", \"Var\", \"Topscoring\", \"Network\", \"All\", \"Composite_Annotated\"). ")
  if (variablespace %in% c("x", "y", "xy") == FALSE)
    stop("Select an adequate variable representation space. x, y, or xy ")
  if (useannot == T && is.null(da_list$annotations)) {
    texts_wizard("\nAnnotation matrix not provided. Considering useannot=F\n")
    texts_wizard( "Please check/repeat da_analysis() including an annotation matrix\n")
    texts_wizard("Read help or follow the tutorial for more information.\n")
    useannot = FALSE
  }
  
    #Reset graphical device  
  if (grDevices::dev.cur() > 1) {
      grDevices::dev.off()
    }
  
    # Required objects ----
  compZ <- NULL
  da_object <- da_list$da
  treatmentstable <- da_list$treatments
  annotations <- da_list$annotations
  datasetnames <- da_list$datasetnames
  treatmentnames <- colnames(treatmentstable)
  originaldata <- da_list$originaldata
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
  
  
  # Top Scoring Loadings plot ----
  if (plottype == "Topscoring") {
    if (useannot == F) {
      topsco_plot <-
        topscoring(
          da_list,
          originaldata = originaldata,
          treatmentstable = treatmentstable,
          n = fortopscoring[2],
          mode = fortopscoring[3],
          numberofcomponents = fortopscoring[1],
          treatment = treatment,
          annots = NULL
        )
    }
    if (useannot == T) {
      topsco_plot <-
        topscoring(
          da_list,
          originaldata = originaldata,
          treatmentstable = treatmentstable,
          n = fortopscoring[2],
          mode = fortopscoring[3],
          numberofcomponents = fortopscoring[1],
          treatment = treatment,
          annots = annotations
        )
    }
    return(topsco_plot)
  }
  
  
  # Distance plot ----
  if (plottype %in% c("Distance","Distance_Ellipse")) {
    if (variablespace == "x") {
      explainedvariance <- as.vector(da_object$explained_variance$X)
      distanceplotdata <- da_object$variates$X
      xaxislabel <-
        paste("Component 1-",
              round(explainedvariance[1] * 100, 2),
              "% of variance",
              sep = " ")
      yaxislabel <-
        paste("Component 2-",
              round(explainedvariance[2] * 100, 2),
              "% of variance",
              sep = " ")
    }
    if (variablespace == "y") {
      explainedvariance <- da_object$explained_variance$Y
      distanceplotdata <- da_object$variates$Y
      xaxislabel <-
        paste("Component 1-",
              round(explainedvariance[1] * 100, 2),
              "% of variance",
              sep = " ")
      yaxislabel <-
        paste("Component 2-",
              round(explainedvariance[2] * 100, 2),
              "% of variance",
              sep = " ")
    }
    if (variablespace == "xy") {
      distanceplotdata <- (da_object$variates$X + da_object$variates$Y) / 2
      xaxislabel <- "XY-Variate 1"
      yaxislabel <- "XY-Variate 2"
    }
    
    if (treatment != 4) {
      # Single factor visualization
      if (treatment == 1)
        tratamientocolores <- as.vector(treatmentstable[, 1])
      if (treatment == 2)
        tratamientocolores <- as.vector(treatmentstable[, 2])
      if (treatment == 3)
        tratamientocolores <- as.vector(treatmentstable[, 3])
      
      tratamientocolores <-
        factor(tratamientocolores, levels = unique(tratamientocolores))
      if(plottype == "Distance_Ellipse"){
        replicates <- nrow(treatmentstable)/length(unique(tratamientocolores))
        if (replicates < 4){plottype = "Distance"}
        if (!is.null(compZ)){plottype = "Distance"}
        
      }
      if(plottype == "Distance_Ellipse"){
        colcolcol <- paleta_tratamientos[1:length(unique(tratamientocolores))]
        if (is.null(compZ)) {# Plot 2D, 1 factor.
          
          h <- ggplot2::ggplot(as.data.frame(distanceplotdata), 
                               ggplot2::aes(distanceplotdata[, 1],
                                            distanceplotdata[, 2],
                                            color = factor(tratamientocolores,levels=unique(tratamientocolores)))) + 
            ggplot2::theme_minimal() + 
            ggplot2::stat_ellipse(level=confidence, 
                                  geom = "polygon", 
                                  alpha=0.3, 
                                  ggplot2::aes(fill=factor(tratamientocolores,levels=unique(tratamientocolores)))) + 
            ggplot2::geom_point()+    
            ggplot2::theme(legend.title = element_blank())
          
          j <- h + ggplot2::scale_color_manual(values = colcolcol) + 
            ggplot2::scale_fill_manual(values=colcolcol) 
          
          da_score_plot3 <- plotly::ggplotly(j)
          
          da_score_plot2 <- plotly::layout(
            da_score_plot3,
            
            xaxis = list(
              title = xaxislabel,
              titlefont = list(size = fontsizes[1]),
              tickfont = list(size = fontsizes[2])
            ),
            yaxis = list(
              title = yaxislabel,
              titlefont = list(size = fontsizes[1]),
              tickfont = list(size = fontsizes[2])
            ),
            title = list(text = "DA Maximum Distance Plot", font = list(size =
                                                                          fontsizes[3])),
            legend = list(font = list(size =
                                        fontsizes[4]))
          )
          
          #retocamos trazas elipses
          da_score_plot2_t <- da_score_plot2
          ellips_text <- paste(confidence*100,"% confidence","<br>interval; ",unique(tratamientocolores))
          
          for(i in (1:length(unique(colcolcol)))){
            da_score_plot2_t <- plotly::style(da_score_plot2_t, 
                                              hoverinfo = "text", text = ellips_text[i], traces = i,name = ellips_text[i], showlegend=T)
          }
          #retocamos trazas samples
          da_score_plot2_tp <- da_score_plot2_t
          track <- c((1+length(unique(colcolcol))):(length(unique(colcolcol))*2)) 
          ptra <- as.data.frame(cbind(rownames(distanceplotdata),treatmentstable[,treatment]))
          a <- lapply(split(ptra, ptra[,2]), function(x) (x[,1]))
          v <- c()
          for(i in 1:length(a)){
            v[i] <- length(a[[i]])
          }
          if(var(v) != 0){
            lenmx <- max(v)
            d <- c()
            for(i in 1:length(a)){
              d[i] <- lenmx - length(a[[i]])
              a[[i]] <- c(a[[i]],rep(NA,d[i]))
            }
          }
          b <- do.call(cbind, a)
          
          for(i in (1:length(unique(colcolcol)))){
            da_distance_plot <- plotly::style(da_score_plot2_tp, 
                                              hoverinfo = "text",
                                              text = paste("Sample: ",stats::na.omit(b[,i]), sep=""),
                                              showlegend = TRUE, 
                                              name = colnames(b)[i],
                                              traces = track[i] )}
          return(da_distance_plot)}
      }else{
      
      da_distance_plot <-
        plotly::plot_ly(
          as.data.frame(distanceplotdata),
          x = ~ distanceplotdata[, 1],
          y = ~ distanceplotdata[, 2],
          hoverinfo = "text",
          text = paste("Sample:", rownames(distanceplotdata), sep =
                         " "),
          showlegend = TRUE,
          type = "scatter",
          mode = "markers",
          color = tratamientocolores,
          colors = paleta_tratamientos[1:length(unique(tratamientocolores))],
          #symbol = factor(rep(c(1:nrow(mcoin_object$mcoa$lambda)),each=nrow(mcoin_object$mcoa$SynVar)), levels=datasetnames),
          #symbols = simbolos[1:length(unique(tratamientocolores))],
          marker = list(size = 11)
        )
      da_distance_plot <-
        plotly::layout(
          p = da_distance_plot,
          xaxis = list(
            title = xaxislabel,
            titlefont = list(size = fontsizes[1]),
            tickfont = list(size = fontsizes[2])
          ),
          yaxis = list(
            title = yaxislabel,
            titlefont = list(size = fontsizes[1]),
            tickfont = list(size = fontsizes[2])
          ),
          title = list(text = "DA Maximum Distance Plot", font = list(size =
                                                                        fontsizes[3])),
          legend = list(font = list(size =
                                      fontsizes[4]))
        )
      
      return(da_distance_plot)}
    }
    
    if (treatment == 4) {
      # 2 factor visualization
      
      da_distance_plot <-
        plotly::plot_ly(
          as.data.frame(distanceplotdata),
          x = ~ distanceplotdata[, 1],
          y = ~ distanceplotdata[, 2],
          hoverinfo = "text",
          text = paste("Sample:", rownames(distanceplotdata), sep =
                         " "),
          showlegend = TRUE,
          type = "scatter",
          mode = "markers",
          color = factor(treatmentstable[, 1], levels = unique(treatmentstable[, 1])),
          colors = paleta_tratamientos[1:length(unique(treatmentstable[, 1]))],
          symbol = factor(treatmentstable[, 2], levels =
                            unique(treatmentstable[, 2])),
          symbols = simbolos[1:length(unique(treatmentstable[, 2]))],
          marker = list(size = 11)
        )
      da_distance_plot <-
        plotly::layout(
          p = da_distance_plot,
          xaxis = list(
            title = xaxislabel,
            titlefont = list(size = fontsizes[1]),
            tickfont = list(size = fontsizes[2])
          ),
          yaxis = list(
            title = yaxislabel,
            titlefont = list(size = fontsizes[1]),
            tickfont = list(size = fontsizes[2])
          ),
          title = list(text = "DA Maximum Distance Plot", font = list(size =
                                                                        fontsizes[3])),
          legend = list(font = list(size =
                                      fontsizes[4]))
        )
      return(da_distance_plot)
      
      
      
    }
  }
  
  # Var plot ----
  if (plottype == "Var") {
    vartable <-
      invisible(mixOmics::plotVar(da_object, plot = F)) # fully functional but to be continued...
    vartable <- vartable[, 1:3]
    if (variablespace == "x") {
      explainedvariance <- as.vector(da_object$explained_variance$X)
      xaxislabel <-
        paste("Component 1-",
              round(explainedvariance[1] * 100, 2),
              "% of variance",
              sep = " ")
      yaxislabel <-
        paste("Component 2-",
              round(explainedvariance[2] * 100, 2),
              "% of variance",
              sep = " ")
    }
    if (variablespace == "y") {
      explainedvariance <- da_object$explained_variance$Y
      xaxislabel <-
        paste("Component 1-",
              round(explainedvariance[1] * 100, 2),
              "% of variance",
              sep = " ")
      yaxislabel <-
        paste("Component 2-",
              round(explainedvariance[2] * 100, 2),
              "% of variance",
              sep = " ")
    }
    if (variablespace == "xy") {
      xaxislabel <- "XY-Variate 1"
      yaxislabel <- "XY-Variate 2"
    }
    
    if (useannot == FALSE) {
      da_VarPlot <- plotly::plot_ly(
        as.data.frame(vartable),
        x = ~ vartable[, 1],
        y = ~ vartable[, 2],
        text = paste("ID:", rownames(vartable)),
        type = "scatter",
        mode = "markers",
        marker = list(size = 6, opacity = 1),
        colors = paleta_niveles[1:length(datasetnames)]
      )
      
      da_VarPlot <-
        plotly::layout(
          da_VarPlot,
          xaxis = list(
            title = xaxislabel,
            titlefont = list(size = fontsizes[1]),
            tickfont = list(size = fontsizes[2])
          ),
          yaxis = list(
            title = yaxislabel,
            titlefont = list(size = fontsizes[1]),
            tickfont = list(size = fontsizes[2])
          ),
          title = list(text = "DA Maximum Distance Plot", font = list(size =
                                                                        fontsizes[3])),
          legend = list(font = list(size = fontsizes[4]))
        )
      
      
      return(da_VarPlot)
    }
    if (useannot == TRUE) {
      # Annotations 
      annotations2 <- as.data.frame(row.names(vartable))
      colnames(annotations2) <- "IDENTIFIER"
      annotations <-
        dplyr::left_join(annotations2, as.data.frame(annotations), by = "IDENTIFIER")
      annotations <- as.matrix(annotations)
      rm(annotations2)
      #FLAG!!!
      # Mapman levels
      mapmanlevels <-
        as.numeric(unique(annotations[order(as.numeric(annotations[, 3])), 3]))
      mapmannames <-
        annotations[order(as.numeric(annotations[, 3])), ]
      mapmannames <- as.vector(unique(mapmannames[, 4]))
      
      da_VarPlot <- plotly::plot_ly(
        as.data.frame(vartable),
        x = ~ vartable[, 1],
        y = ~ vartable[, 2],
        text = paste(
          "ID:",
          rownames(vartable),
          "</br> Description: ",
          annotations[, 2],
          "</br> Mercator Bin: ",
          annotations[, 4]
        ),
        type = "scatter",
        mode = "markers",
        marker = list(size = 6, opacity = 1),
        color = factor(annotations[, 4], levels =
                         mapmannames),
        colors = paleta_mapman[mapmanlevels]
      )
      
      da_VarPlot <-
        plotly::layout(
          da_VarPlot,
          xaxis = list(
            title = xaxislabel,
            titlefont = list(size = fontsizes[1]),
            tickfont = list(size = fontsizes[2])
          ),
          yaxis = list(
            title = yaxislabel,
            titlefont = list(size = fontsizes[1]),
            tickfont = list(size = fontsizes[2])
          ),
          title = list(text = "DA Maximum Distance Plot", font = list(size =
                                                                        fontsizes[3])),
          legend = list(font = list(size = fontsizes[4]))
        )
      
      return(da_VarPlot)
    }
  }
  # Composite plot ----
  if (plottype == "Composite") {
    p1 <-
      da_plot(da_list,
              treatment,
              plottype = "Distance",
              fontsizes = fontsizes)
    p2 <-
      da_plot(da_list,
              treatment,
              plottype = "Var",
              fontsizes = fontsizes)
    da_CompositeplotPlot <-
      plotly::subplot(
        list(p1, p2),
        nrows = 2,
        shareX = F,
        shareY = F,
        margin = 0.04
      )
    da_CompositeplotPlot <- plotly::layout(da_CompositeplotPlot,
                                           title = list(text = "DA Analysis", font = list(size =
                                                                                            fontsizes[3])))
    return(da_CompositeplotPlot)
  }
  # Composite Annotated ----
  if (plottype == "Composite_Annotated") {
    p1 <-
      da_plot(da_list,
              treatment,
              plottype = "Distance",
              fontsizes = fontsizes)
    p2 <-
      da_plot(
        da_list,
        treatment,
        plottype = "Var",
        fontsizes = fontsizes,
        useannot = T
      )
    da_CompositeplotPlot <-
      plotly::subplot(
        list(p1, p2),
        nrows = 2,
        shareX = F,
        shareY = F,
        margin = 0.04
      )
    da_CompositeplotPlot <- plotly::layout(
      da_CompositeplotPlot,
      title = list(text = "DA Analysis", font = list(size =
                                                       fontsizes[3])),
      showlegend = FALSE
    )
    return(da_CompositeplotPlot)
  }
  # Network plot ----
  if (plottype == "Network") {
    if (useannot == FALSE) {
      set.seed(5881)
      
      what <- tryCatch(suppressMessages(mixOmics::network(da_object, comp = c(1,2), cutoff = cutoff)), error=function(e) {
        while(dev.cur()>1){dev.off()}
        texts_wizard("\nThe obtained network is too large for be displayed in RStudio plot viewer")
        texts_wizard("\nBuuuuut, it's ok, we have exported network.gml file compatible with cytoscape, and network.pdf in your working directory instead")
        #yeah! we are that good & kind
        #you are welcome
        grDevices::pdf("da_network.pdf")
        suppressMessages(mixOmics::network(da_object, cutoff = cutoff, comp = c(1, 2)))
        dev.off() 
        grDevices::pdf("da_network.pdf")
        a <- suppressMessages(mixOmics::network(da_object,  cutoff = cutoff, comp = c(1, 2)))
        dev.off()
        igraph::write.graph(a$gR, file = "da_network.gml", format = "gml")
        class(da_object) <- c(class(da_object), "networkplot")
        
      }
      
      )
      return(da_object)
    }
    if (useannot == TRUE) {
      
      
      what <- tryCatch(suppressMessages(mixOmics::network(da_object, comp = c(1, 2), cutoff = cutoff)), error=function(e) {
        while(dev.cur()>1){dev.off()}
        texts_wizard("\nThe obtained network is too large for be displayed in RStudio plots")
        texts_wizard("\nBuuuuut, it's ok, we have exported network.gml file compatible with cytoscape, and  network.pdf in your working directory instead")
        #yeah! we are that good & kind
        #you are welcome
        grDevices::pdf("da_network.pdf")
        suppressMessages(mixOmics::network( da_list[[5]],
                           comp = c(1, 2),cex.node.name = 0.5,cutoff = cutoff))          
        dev.off() 
        grDevices::pdf("da_network.pdf")
        a <-suppressMessages(mixOmics::network(da_list[[5]],comp = c(1, 2),cex.node.name = 0.5,cutoff = cutoff))
        dev.off()
        igraph::write.graph(a$gR, file = "da_network.gml", format = "gml")
        class(da_list[[5]]) <-
          c(class(da_list[[5]]), "networkplotannot")
      }
      )
      return(da_list[[5]])
    }
  }
  
}
