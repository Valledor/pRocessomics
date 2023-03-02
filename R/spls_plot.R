

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


#' @name spls_plot
#' @title sparse Partial Least Squares plot
#' @description A function to depict sPLS analysis
#' @usage spls_plot(spls_list, treatment = 1, plottype = "Distance",
#' variablespace = "xy", useannot = FALSE, fontsizes = c(14,10,16,12),
#' cutoff = 0.8, fortopscoring = c(1,20,"abs"), confidence = 0.9)
#' @param spls_list spls analysis object, class "splsanalysis", see pRocessomics::spls_analysis function
#' @param treatment numeric value to determine the colors and symbols to be displayed in the plot, used in Score plot and Topscring plot. Available optiones are:
#' \itemize{
#' \item 1 to show different colors according to treatment 1 (treatment1col)
#' \item 2 to show different colors according to treatment 2 (treatment2col)
#' \item 3 to show different colors according to the combination of both treatments
#' \item 4 to show different colors according to treatment 1 and different symbols for treatment 2
#' }
#' @param plottype available options:
#' \itemize{
#' \item Distance; sample plot each dot represent a sample
#' \item Distance_Ellipse, distance plot with confidence ellipses by treatment
#' \item Topscoring; most representantive loadings according to Distance plot
#' \item CompoX-Y; composition of Distance plots of X and Y variablespace (separately)
#' \item Var; variable plot, each dot represent a variable. select dataset and annotation
#' \item Composite; combination of the Distance (X, Y and XY variablespaces separately) and Var plot
#' \item Network; network
#' \item All; all of above
#' }
#' @param variablespace c("x","y","xy")
#' @param useannot logical, indicating if annotation should be used or not, for annotations use the splsda analysis should include it
#' @param fontsizes font sizes for plot elements
#' @param cutoff cutoff value for network prunning
#' @param fortopscoring vector indicating the number of the component to be analyzed, the top n loadings to plot and the mode to compute "abs" for absolute value and "posneg" for positive and negative loading analyzed separately
#' @param confidence numeric value between 0 and 1 for confidence ellipses in Distance_Ellipse plot
#' @return plotly type plot
#' @author Luis Valledor and Laura Lamelas

#' @export
#' @importFrom plotly plot_ly layout subplot style ggplotly
#' @importFrom mixOmics plotVar network
#' @importFrom dplyr left_join
#' @importFrom igraph write.graph
#' @importFrom grDevices pdf
#' @importFrom ggplot2 ggplot theme scale_color_manual scale_fill_manual stat_ellipse geom_point theme_minimal aes
#' @importFrom stats na.omit



spls_plot <-
  function(spls_list,
           treatment = 1,
           plottype = "Distance",
           variablespace = "xy",
           useannot = FALSE,
           fontsizes = c(14, 10, 16, 12),
           cutoff = 0.8,
           fortopscoring = c(1, 20, "abs"),
           confidence = 0.9) {
    # Checks ----
   
      if(!inherits(spls_list, "splsanalysis"))
      stop("SPLS analysis object is expected, please run spls_analysis first")
    if (plottype %in% c(
      "Composite",
      "CompoX-Y",
      "Distance",
      "Distance_Ellipse",
      "Var",
      "Network",
      "Topscoring",
      "All",
      "Composite_Annotated"
    ) == FALSE)
      stop("Please select a valid plot type")
    if (variablespace %in% c("x", "y", "xy") == FALSE)
      stop("Select an adequate variable representation space. x, y, or xy ")
    if (useannot == T && is.null(spls_list$annotations)) {
      texts_wizard("\nAnnotation matrix not provided. Considering useannot=F\n")
      texts_wizard("Please check/repeat spls_analysis() including an annotation matrix\n")
      texts_wizard("Read help or follow the tutorial for more information.\n")
      useannot = FALSE
    }
    
    #Reset graphical device  
    if (grDevices::dev.cur() > 1) {
      grDevices::dev.off()
    }
    
    # Initial objects ----
    compZ <- NULL
    spls_object <- spls_list$splso
    treatmentstable <- spls_list$treatments
    annotations <- spls_list$annotations
    datasetnames <- spls_list$datasetnames
    treatmentnames <- colnames(treatmentstable)
    originaldata <- spls_list$originaldata
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
    
    
    # Topscoring plot ----
    if (plottype == "Topscoring") {
      if (useannot == T) {
        topsco_plot <-
          topscoring(
            spls_list,
            originaldata = originaldata,
            treatmentstable = treatmentstable,
            n = fortopscoring[2],
            mode = fortopscoring[3],
            numberofcomponents = fortopscoring[1],
            treatment = treatment,
            annots = annotations
          )
      }
      if (useannot == F) {
        topsco_plot <-
          topscoring(
            spls_list,
            originaldata = originaldata,
            treatmentstable = treatmentstable,
            n = fortopscoring[2],
            mode = fortopscoring[3],
            numberofcomponents = fortopscoring[1],
            treatment = treatment,
            annots = NULL
          )
      }
      return(topsco_plot)
    }
    # Distance plot ----
    if (plottype %in% c("Distance","Distance_Ellipse")) {
      if (variablespace == "x") {
        explainedvariance <- as.vector(spls_object$explained_variance$X)
        distanceplotdata <- spls_object$variates$X
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
        explainedvariance <- spls_object$explained_variance$Y
        distanceplotdata <- spls_object$variates$Y
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
        distanceplotdata <-
          (spls_object$variates$X + spls_object$variates$Y) / 2
        xaxislabel <- "XY-Variate 1"
        yaxislabel <- "XY-Variate 2"
      }
      
      if (treatment != 4) {
        #Asi seleccionamos todas las visualizaciones que solo tengan un factor a mostrar.
        if (treatment == 1)
          tratamientocolores <- as.vector(treatmentstable[, 1])
        if (treatment == 2)
          tratamientocolores <- as.vector(treatmentstable[, 2])
        if (treatment == 3)
          tratamientocolores <- as.vector(treatmentstable[, 3])
        
        tratamientocolores <- factor(tratamientocolores, levels = unique(tratamientocolores))
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
            
            spls_score_plot3 <- plotly::ggplotly(j)
            
            spls_score_plot2 <- plotly::layout(
              spls_score_plot3,
              
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
              title = list(text = "SPLS Maximum Distance Plot", font = list(size =
                                                                            fontsizes[3])),
              legend = list(font = list(size =
                                          fontsizes[4]))
            )
            
            #retocamos trazas elipses
            spls_score_plot2_t <- spls_score_plot2
            ellips_text <- paste(confidence*100,"% confidence","<br>interval; ",unique(tratamientocolores))
            
            for(i in (1:length(unique(colcolcol)))){
              spls_score_plot2_t <- plotly::style(spls_score_plot2_t, 
                                                  hoverinfo = "text", text = ellips_text[i], traces = i,name = ellips_text[i], showlegend=T)
            }
            #retocamos trazas samples
            spls_score_plot2_tp <- spls_score_plot2_t
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
              spls_distance_plot <- plotly::style(spls_score_plot2_tp, 
                                                  hoverinfo = "text",
                                                  text = paste("Sample: ",stats::na.omit(b[,i]), sep=""),
                                                  showlegend = TRUE, 
                                                  name = colnames(b)[i],
                                                  traces = track[i] )}
            return(spls_distance_plot)}
        }else{
         spls_distance_plot <-
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
        spls_distance_plot <-
          plotly::layout(
            p = spls_distance_plot,
            xaxis = list(
              title = list(
                text = xaxislabel,
                font = list(size = fontsizes[1])
                ),
              tickfont = list(size = fontsizes[2])
            ),
            yaxis = list(
              title = list(
                text = yaxislabel,
                font = list(size = fontsizes[1])
                ),
              tickfont = list(size = fontsizes[2])
            ),
            title = list(
              text = "SPLS Maximum Distance Plot",
              font = list(size = fontsizes[3])
              ),
            legend = list(font = list(size = fontsizes[4]))
          )
        
        return(spls_distance_plot)}
      }
      
      if (treatment == 4) {
        #Asi seleccionamos las visualizaciones con 2 factores a mostrar
        
        spls_distance_plot <-
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
        spls_distance_plot <-
          plotly::layout(
            p = spls_distance_plot,
            xaxis = list(
              title = list(
                text = xaxislabel,
                font = list(size = fontsizes[1])
                ),
              tickfont = list(size = fontsizes[2])
            ),
            yaxis = list(
              title = list(
                text = yaxislabel,
                font = list(size = fontsizes[1])
                ),
              tickfont = list(size = fontsizes[2])
            ),
            title = list(
              text = "SPLS Maximum Distance Plot",
              font = list(size = fontsizes[3])
              ),
            legend = list(font = list(size = fontsizes[4]))
          )
        return(spls_distance_plot)
      }
    }
    # Var plot ----
    if (plottype == "Var") {
      vartable <-
        invisible(mixOmics::plotVar(spls_object, plot = FALSE)) #to be replaced
      vartable <- vartable[, 1:3]
      if (variablespace == "x") {
        explainedvariance <- as.vector(spls_object$explained_variance$X)
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
        explainedvariance <- spls_object$explained_variance$Y
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
        spls_VarPlot <- plotly::plot_ly(
          as.data.frame(vartable),
          x = ~ vartable[, 1],
          y = ~ vartable[, 2],
          text = paste("ID:", rownames(vartable)),
          type = "scatter",
          mode = "markers",
          marker = list(size = 6, opacity = 1),
          colors = paleta_niveles[1:length(datasetnames)],
          symbol = factor(vartable[, 3], levels = unique(vartable[, 3])),
          symbols = simbolos[1:length(unique(vartable[, 3]))]
        )
        
        spls_VarPlot <-
          plotly::layout(
            spls_VarPlot,
            xaxis = list(
              title = list(
                text = xaxislabel,
                font = list(size = fontsizes[1])
                ),
              tickfont = list(size = fontsizes[2])
            ),
            yaxis = list(
              title = list(
                text= yaxislabel,
                font = list(size = fontsizes[1])
                ),
              tickfont = list(size = fontsizes[2])
            ),
            title = list(
              text = "sPLS Variable Plot",
              font = list(size = fontsizes[3])
              ),
            legend = list(font = list(size = fontsizes[4]))
          )
        return(spls_VarPlot)
      }
      
      if (useannot == TRUE) {
        #FLAG
        
        annotations2 <- as.data.frame(row.names(vartable))
        colnames(annotations2) <- "IDENTIFIER"
        annotations <-
          dplyr::left_join(annotations2, as.data.frame(annotations), by = "IDENTIFIER") #con esto asignamos cada anotación a la variable correspondiente por rowname. Asi no hay fallo.
        annotations <-
          as.matrix(annotations) #Por alguna razon, sin esto, luego no se ordenan los factores ni los colores.
        rm(annotations2)
        
        #Sacamos los niveles del mapman, y las etiquetas
        mapmanlevels <-
          as.numeric(unique(annotations[order(as.numeric(annotations[, 3])), 3]))
        mapmannames <-
          annotations[order(as.numeric(annotations[, 3])), ]
        mapmannames <- as.vector(unique(mapmannames[, 4]))
        
        spls_VarPlot <- plotly::plot_ly(
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
          colors = paleta_mapman[mapmanlevels],
          symbol = factor(vartable[, 3], levels = unique(vartable[, 3])),
          symbols = simbolos[1:length(unique(vartable[, 3]))]
        )
        
        spls_VarPlot <-
          plotly::layout(
            spls_VarPlot,
            xaxis = list(
              title = list(
                text = xaxislabel,
                font = list(size = fontsizes[1])
                ),
              tickfont = list(size = fontsizes[2])
            ),
            yaxis = list(
              title = list(
                text = yaxislabel,
                font = list(size = fontsizes[1])
                ),
              tickfont = list(size = fontsizes[2])
            ),
            title = list(
              text = "sPLS Variable Plot",
              font = list(size = fontsizes[3])
              ),
            legend = list(font = list(size = fontsizes[4]))
          )
        
        return(spls_VarPlot)
      }
    }
    # Compo X-Y plot ----
    if (plottype == "CompoX-Y") {
      p1 <-
        spls_plot(
          spls_list,
          treatment,
          plottype = "Distance",
          fontsizes = fontsizes,
          variablespace = "x"
        )
      p2 <-
        spls_plot(
          spls_list,
          treatment,
          plottype = "Distance",
          fontsizes = fontsizes,
          variablespace = "y"
        )
      spls_CompositeplotPlot <-
        plotly::subplot(
          plotly::subplot(p1, p2),
          nrows = 2,
          shareX = F,
          shareY = F,
          margin = 0.04
        )
      spls_CompositeplotPlot <-
        plotly::layout(spls_CompositeplotPlot,
                       title = list(
                         text = "SPLS Analysis", 
                         font = list(size = fontsizes[3])
                         )
                       )
      return(spls_CompositeplotPlot)
    }
    # Composite plot ----
    if (plottype == "Composite") {
      p1 <-
        spls_plot(
          spls_list,
          treatment,
          plottype = "Distance",
          fontsizes = fontsizes,
          variablespace = "x"
        )
      p2 <-
        spls_plot(
          spls_list,
          treatment,
          plottype = "Distance",
          fontsizes = fontsizes,
          variablespace = "y"
        )
      p3 <-
        spls_plot(
          spls_list,
          treatment,
          plottype = "Distance",
          fontsizes = fontsizes,
          variablespace = "xy"
        )
      p4 <-
        spls_plot(spls_list,
                  treatment,
                  plottype = "Var",
                  fontsizes = fontsizes)
      spls_CompositeplotPlot <-
        plotly::subplot(
          plotly::subplot(p1, p2),
          plotly::subplot(p3, p4),
          nrows = 2,
          shareX = F,
          shareY = F,
          margin = 0.04
        )
      spls_CompositeplotPlot <-
        plotly::layout(spls_CompositeplotPlot,
                       title = list(
                         text = "SPLS Analysis", 
                         font = list(size = fontsizes[3])
                         )
                       )
      return(spls_CompositeplotPlot)
    }
    # Composite_Annotated plot ----
    if (plottype == "Composite_Annotated") {
      p1 <-
        spls_plot(
          spls_list,
          treatment,
          plottype = "Distance",
          fontsizes = fontsizes,
          variablespace = "x"
        )
      p2 <-
        spls_plot(
          spls_list,
          treatment,
          plottype = "Distance",
          fontsizes = fontsizes,
          variablespace = "y"
        )
      p3 <-
        spls_plot(
          spls_list,
          treatment,
          plottype = "Distance",
          fontsizes = fontsizes,
          variablespace = "xy"
        )
      p4 <-
        spls_plot(
          spls_list,
          treatment,
          plottype = "Var",
          fontsizes = fontsizes,
          useannot = T
        )
      spls_CompositeplotPlot <-
        plotly::subplot(
          plotly::subplot(p1, p2),
          plotly::subplot(p3, p4),
          nrows = 2,
          shareX = F,
          shareY = F,
          margin = 0.04
        )
      spls_CompositeplotPlot <-
        plotly::layout(
          p = spls_CompositeplotPlot,
          title = list(text = "SPLS Analysis", font = list(size = fontsizes[3])),
          showlegend = FALSE
        )
      return(spls_CompositeplotPlot)
    }
    # Network plot ----
    if (plottype == "Network") {
      if (useannot == FALSE) {
        set.seed(5881)
        
        what <- tryCatch(suppressMessages(mixOmics::network(spls_object, comp = c(1, 2), cutoff = cutoff)), error=function(e) {
          while(dev.cur()>1){dev.off()}
          texts_wizard("\nThe obtained network is too large for be displayed in RStudio plot viewer")
          texts_wizard("\nBuuuuut, it's ok, we have exported network.gml file compatible with cytoscape, and  network.pdf in your working directory instead")
          #yeah! we are that good & kind
          #you are welcome
          grDevices::pdf("network.pdf")
          suppressMessages(mixOmics::network(spls_object, comp = c(1, 2), cutoff = cutoff))
          dev.off() 
          grDevices::pdf("network.pdf")
          a <- suppressMessages(mixOmics::network(spls_object, comp = c(1, 2), cutoff = cutoff))
          dev.off()
          igraph::write.graph(a$gR, file = "network.gml", format = "gml")
          class(spls_object) <- c(class(spls_object), "networkplot")
          
        }
        
        )
        
          
       
        
        
        
        
        
        return(spls_object)
      }
      if (useannot == TRUE) {
        
        
        what <- tryCatch(suppressMessages(mixOmics::network(spls_list[[5]], comp = c(1, 2), cutoff = cutoff)), error=function(e) {
          while(dev.cur()>1){dev.off()}
          texts_wizard("\nThe obtained network is too large for be displayed in RStudio plot viewer")
          texts_wizard("\nBuuuuut, it's ok, we have exported network.gml file compatible with cytoscape, and  network.pdf in your working directory instead")
          #yeah! we are that good & kind
          #you are welcome
          grDevices::pdf("network.pdf")
          suppressMessages(mixOmics::network(spls_list[[5]],comp = c(1, 2),cex.node.name = 0.5,cutoff = cutoff))        
          dev.off() 
          grDevices::pdf("network.pdf")
          a <-suppressMessages(mixOmics::network(spls_list[[5]],comp = c(1, 2),cex.node.name = 0.5,cutoff = cutoff))
          dev.off()
          igraph::write.graph(a$gR, file = "network.gml", format = "gml")
          class(spls_list[[5]]) <-c(class(spls_list[[5]]), "networkplotannot")
        }
        )
      return(spls_list[[5]])
      }
    }
    # All ----
  }
