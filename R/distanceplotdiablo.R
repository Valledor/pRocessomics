#' @importFrom plotly plot_ly layout style ggplotly
#' @importFrom ggplot2 ggplot theme scale_color_manual scale_fill_manual stat_ellipse geom_point theme_minimal aes
#' @importFrom stats na.omit


distanceplotdiablo <-
  function(variates,
           varianceexplained,
           treatment,
           treatmentstable,
           levelname,
           fontsizes = c(14, 10, 16, 12),
           confidence,
           plottype = "Distance_Ellipse") {
    
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
    
    
    explainedvariance <- as.vector(varianceexplained)
    
    compZ <- NULL
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
    
    if (treatment != 4) {
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
          
          h <- ggplot2::ggplot(as.data.frame(variates), 
                               ggplot2::aes(variates[, 1],
                                            variates[, 2], 
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
          
          bda_score_plot3 <- plotly::ggplotly(j)
          
          bda_score_plot2 <- plotly::layout(
            bda_score_plot3,
            
            xaxis = list(
              title = list(text = xaxislabel,
                           font = list(size = fontsizes[1])
              ),
              tickfont = list(size = fontsizes[2])
            ),
            yaxis = list(
              title = list (text = yaxislabel,
                            font = list(size = fontsizes[1])),
              tickfont = list(size = fontsizes[2])
            ),
            annotations = list( x = 0, 
                                y = 1.02, 
                                text = levelname, 
                                font = list(size = fontsizes[3]),
                                showarrow = F,
                                xref = 'paper',
                                yref = 'paper'
            )
          )
          
          #retocamos trazas elipses
          bda_score_plot2_t <- bda_score_plot2
          ellips_text <- paste(confidence*100,"% confidence","<br>interval; ",unique(tratamientocolores))
          
          for(i in (1:length(unique(colcolcol)))){
            bda_score_plot2_t <- plotly::style(bda_score_plot2_t, 
                                               hoverinfo = "text", text = ellips_text[i], traces = i,name = ellips_text[i], showlegend=T)
          }
          #retocamos trazas samples
          bda_score_plot2_tp <- bda_score_plot2_t
          track <- c((1+length(unique(colcolcol))):(length(unique(colcolcol))*2)) 
          ptra <- as.data.frame(cbind(rownames(variates),treatmentstable[,treatment]))
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
            da_distance_plot <- plotly::style(bda_score_plot2_tp, 
                                                hoverinfo = "text",
                                                text = paste("Sample: ",stats::na.omit(b[,i]), sep=""),
                                                showlegend = TRUE, 
                                                name = colnames(b)[i],
                                                traces = track[i] )}
          return(da_distance_plot)}
      }
      
      da_distance_plot <- plotly::plot_ly(
        as.data.frame(variates),
        x = ~ variates[, 1],
        y = ~ variates[, 2],
        hoverinfo = "text",
        text = paste("Sample:", rownames(variates), sep =
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
          da_distance_plot,
          xaxis = list(
            title = list(text = xaxislabel,
                        font = list(size = fontsizes[1])
            ),
            tickfont = list(size = fontsizes[2])
          ),
          yaxis = list(
            title = list (text = yaxislabel,
            font = list(size = fontsizes[1])),
            tickfont = list(size = fontsizes[2])
          ),
          annotations = list( x = 0, 
                              y = 1.02, 
                              text = levelname, 
                              font = list(size = fontsizes[3]),
                              showarrow = F,
                              xref = 'paper',
                              yref = 'paper'
                              )
          #legend = list(font = list(size = fontsizes[4]))
        )
      
      return(da_distance_plot)
    }
    
    if (treatment == 4) {
      #Asi seleccionamos las visualizaciones con 2 factores a mostrar
      
      da_distance_plot <- plotly::plot_ly(
        as.data.frame(variates),
        x = ~ variates[, 1],
        y = ~ variates[, 2],
        hoverinfo = "text",
        text = paste("Sample:", rownames(variates), sep =
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
          #title = list(text = levelname,
          #             font = list(size = fontsizes[3])),
          #legend = list(font = list(size = fontsizes[4]))
          annotations = list( x = 0.0, 
                              y = 0.9, 
                              text = levelname, 
                              font = list(size = fontsizes[3]),
                              showarrow = F,
                              xref = 'paper',
                              yref = 'paper'
          )
        )
      return(da_distance_plot)
      
    }
  }
