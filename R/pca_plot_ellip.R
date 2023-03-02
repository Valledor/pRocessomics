

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


#' @name pca_plot
#' @title PCA ploting function
#' @description A function to plot Principal Component Analysis (PCA) results
#' @usage pca_plot(pca_list,treatment = 1, compX = 1, compY = 2, compZ = NULL, 
#' plottype = "Composite", useannot = FALSE, fontsizes = c(14,10,16,12), 
#' fortopscoring = c(1,20,"abs"), confidence = 0.9)
#' @param pca_list PCA object to plot, created with pRocessomics::pca_analysis() function
#' @param treatment Numeric value to determine the colors and symbols to be displayed in the plot, used in Score or Score_ellipse plot, Topscring plot and Biplot. Available options are:
#' \itemize{
#' \item 1 to show different colors according to treatment 1 (treatment1col)
#' \item 2 to show different colors according to treatment 2 (treatment2col)
#' \item 3 to show different colors according to the combination of both treatments
#' \item 4 to show different colors according to treatment 1 and different symbols for treatment 2
#' }
#' @param compX component to be plotted in X axis in Score and Var plots, default 1
#' @param compY component to be plotted in X axis in Score and Var plots, default 2
#' @param compZ component to be plotted in Z axis for 3D plots in Score and Var plots, default NULL
#' @param plottype available options:
#' \itemize{
#' \item Scree, explained variance plot
#' \item Score, sample plot each dot represent a sample
#' \item Score_Ellipse, score plot with confidence ellipses by treatment
#' \item Topscoring, top scoring variables in the selected (PC) principal component
#' \item Biplot, PCA biplot showing samples and variables
#' \item Var, variable plot, each dot represent a variable
#' \item Composite: combination of the most representative plots (Score, Var and Biplot)
#' \item All all aforementioned plots will be exported as .pdf files in your working directory
#' }
#' @param useannot logical, indicating if annotations should be used or not, please note the previous pca analysis should include annotations
#' @param confidence numeric value between 0 and 1 for confidence ellipses in Score_Ellipse plot
#' @param fortopscoring vector indicating the number of the component to be analyzed, the top n loadings to plot and the mode to compute "abs" for absolute value and "posneg" for positive and negative loading analyzed separately
#' @param fontsizes vector containing the font sizes to use in the plots
#' @return pca plot
#' @author Luis Valledor and Laura Lamelas

#' @export
#' @importFrom plotly plot_ly subplot layout add_trace add_markers style ggplotly
#' @importFrom ggplot2 ggplot theme scale_color_manual scale_fill_manual stat_ellipse geom_point theme_minimal aes
#' @importFrom dplyr left_join
#' @importFrom stats na.omit

pca_plot <-
  function(pca_list,
           treatment = 1,
           compX = 1,
           compY = 2,
           compZ = NULL,
           plottype = "Composite",
           useannot = FALSE,
           fontsizes = c(14, 10, 16, 12),
           fortopscoring = c(1, 20, "abs"),
           confidence = 0.9) {
    # Checks ----
    options(warn = -1)
  
    if (!inherits(pca_list, "pcaanalysis"))
      stop("PCA analysis object is expected, please run pca_analysis first")
    if (!inherits(pca_list[[1]], "prcomp"))
      stop("PCA analysis object seems to be wrong. Please re-run pca_analysis")
    if (plottype %in% c("Composite", "Score", "Scree", "Var", "Biplot", "Topscoring", "All", "Score_Ellipse") == FALSE)
      stop(
        'Please select a valid plot type ("Composite","Score","Score_Ellipse","Scree","Var","Biplot","Topscoring","All")'
      )
    if (useannot == T && is.null(pca_list$annotations)) {
      cat(
        paste(
          "\nAnnotation matrix not provided. Considering useannot=F\n",
          "Please check/repeat pca_analysis() including an annotation matrix\n",
          "Read help or follow the tutorial for more information.\n",
          sep = ""
        )
      )
      useannot = FALSE
    }
    
    # Required objects ----
    
    pca_object <- pca_list$pca
    treatmentstable <- pca_list$treatments
    annotations <- pca_list$annotations
    datasetnames <- pca_list$datasetnames
    treatmentnames <- colnames(treatmentstable)
    originaldata <- pca_list$originaldata
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
    
    
    # Explained variance ----
    explainedvariance <-
      round((pca_object$sdev ^ 2) * 100 / sum(pca_object$sdev ^ 2), 2)  #Expalined variance calculation by each component
    
    # Topscoring plot ----
    if (plottype == "Topscoring") {
      if (useannot == F) {
        topsco_plot <-
          topscoring(
            pca_list,
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
            pca_list,
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
    
    # Score plot or Score_ellipse----
    if (plottype %in% c("Score","Score_Ellipse")) {
      
      if (compX > ncol(pca_object$x))
        stop(
          paste(
            "\nX-Axis component do not exist. Only ",
            ncol(pca_object$x),
            " components were defined in PCA analysis",
            sep = ""
          )
        )
      if (compY > ncol(pca_object$x))
        stop(
          paste(
            "\nY-Axis component do not exist. Only ",
            ncol(pca_object$x),
            " components were defined in PCA analysis",
            sep = ""
          )
        )
      
      if (treatment != 4) {
        if (treatment == 1)
          tratamientocolores <- as.vector(treatmentstable[, 1])
        if (treatment == 2)
          tratamientocolores <- as.vector(treatmentstable[, 2])
        if (treatment == 3)
          tratamientocolores <- as.vector(treatmentstable[, 3])
        
        tratamientocolores <-
          factor(tratamientocolores, levels = unique(tratamientocolores))
        
        if(plottype == "Score_Ellipse"){
          replicates <- nrow(treatmentstable)/length(unique(tratamientocolores))
          if (replicates < 4){plottype = "Score"}
          if (!is.null(compZ)){plottype = "Score"}
         
        }
        if(plottype == "Score_Ellipse"){
          colcolcol <- paleta_tratamientos[1:length(unique(tratamientocolores))]
        if (is.null(compZ)) {# Plot 2D, 1 factor.
         
          h <- ggplot2::ggplot(as.data.frame(pca_object$x), 
                               ggplot2::aes(pca_object$x[, compX],
                                   pca_object$x[, compY], 
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
                 
          pca_score_plot3 <- plotly::ggplotly(j)
          
          pca_score_plot2 <- plotly::layout(
            pca_score_plot3,
            
            xaxis = list(
              zeroline=T,
              title = list(
                text = paste(
                  "Component",
                  compX,
                  "-",
                  explainedvariance[compX],
                  "% of variance",
                  sep = " "
                ),
                font =list(
                  size = fontsizes[1]
                  )
              ),
              tickfont = list(size =
                                fontsizes[2])
            ),
            yaxis = list(
              zeroline=T,
              title = list(
                text = paste(
                  "Component",
                  compY,
                  "-",
                  explainedvariance[compY],
                  "% of variance",
                  sep = " "
                ),
                font =
                  list(size = fontsizes[1])
              ),
              tickfont = list(size =
                                fontsizes[2])
            ),
            title = list(text = "PCA Score Plot", font =
                           list(size = fontsizes[3])),
            legend = list(font = list(size =
                                         fontsizes[4]))
          )
          
          #retocamos trazas elipses
          pca_score_plot2_t <- pca_score_plot2
          ellips_text <- paste(confidence*100,"% confidence","<br>interval; ",unique(tratamientocolores))
         
          for(i in (1:length(unique(colcolcol)))){
          pca_score_plot2_t <- plotly::style(pca_score_plot2_t, 
                                          hoverinfo = "text", text = ellips_text[i], traces = i,name = ellips_text[i], showlegend=T)
                                            }
          #retocamos trazas samples
          pca_score_plot_tp <- pca_score_plot2_t
          track <- c((1+length(unique(colcolcol))):(length(unique(colcolcol))*2)) 
          ptra <- as.data.frame(cbind(rownames(pca_object$x),treatmentstable[,treatment]))
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
          pca_score_plot_tp <- plotly::style(pca_score_plot_tp, 
                                          hoverinfo = "text",
                                          text = paste("Sample: ",stats::na.omit(b[,i]), sep=""),
                                          showlegend = TRUE, 
                                          name = colnames(b)[i],
                                          traces = track[i] )}
            return(pca_score_plot_tp)}
        }
        else{
          pca_score_plot <-
            plotly::plot_ly(
              as.data.frame(pca_object$x),
              x = ~ pca_object$x[, compX],
              y = ~ pca_object$x[, compY],
              hoverinfo = "text",
              text = paste("Sample:", rownames(pca_object$x), sep =" "),
              showlegend = TRUE,
              type = "scatter",
              mode = "markers",
              color = tratamientocolores,
              colors = paleta_tratamientos[1:length(unique(tratamientocolores))],
              #symbol = tratamientocolores,
              #symbols = simbolos[1:length(unique(tratamientocolores))],
              marker = list(size = 11)
            )
          pca_score_plot <- plotly::layout(
            pca_score_plot,
            
            xaxis = list(
              title = list(
                text = paste(
                  "Component",
                  compX,
                  "-",
                  explainedvariance[compX],
                  "% of variance",
                  sep = " "
                ),
                font =list(
                  size = fontsizes[1]
                )
              ),
              tickfont = list(size =
                                fontsizes[2])
            ),
            yaxis = list(
              title = list(
                text = paste(
                  "Component",
                  compY,
                  "-",
                  explainedvariance[compY],
                  "% of variance",
                  sep = " "
                ),
                font =
                  list(size = fontsizes[1])
              ),
              tickfont = list(size =
                                fontsizes[2])
            ),
            title = list(text = "PCA Score Plot", font =
                           list(size = fontsizes[3])),
            legend = list(font = list(size =
                                        fontsizes[4]))
          )
          return(pca_score_plot)
        }
        
        if (!is.null(compZ)) {
          # Plot 3D, 1 factor.
          if(plottype == "Score_ellipse") print("\nSorry! Confidence ellipses are not implemented for 3D plots yet")
          if (compZ > ncol(pca_object$x))
            stop(
              paste(
                "\nZ-Axis component do not exist. Only ",
                ncol(pca_object$x),
                " components were defined in PCA analysis",
                sep = ""
              )
            )
          pca_score_plot <-
            plotly::plot_ly(
              as.data.frame(pca_object$x),
              x = ~ pca_object$x[, compX],
              y = ~ pca_object$x[, compY],
              z = ~ pca_object$x[, compZ],
              hoverinfo = "text",
              text = paste("Sample:", rownames(pca_object$x), sep =
                             " "),
              showlegend = TRUE,
              type = "scatter3d",
              mode = "markers",
              color = tratamientocolores,
              colors = paleta_tratamientos[1:length(unique(tratamientocolores))],
              #symbol = tratamientocolores,
              #symbols = simbolos[1:length(unique(tratamientocolores))],
              marker = list(size = 11)
            )
          pca_score_plot <- plotly::layout(
            pca_score_plot,
            scene = list(
              xaxis = list(
                title = list(
                  text = (
                    paste(
                      "Component",
                      compX,
                      "-",
                      explainedvariance[compX],
                      "% of variance",
                      sep = " "
                    )
                  ),
                  font = list(size = fontsizes[1])
                ),
                tickfont = list(size = fontsizes[2])
              ),
              yaxis = list(
                title = list(
                  text = (
                    paste(
                      "Component",
                      compY,
                      "-",
                      explainedvariance[compY],
                      "% of variance",
                      sep = " "
                    )
                  ),
                  font = list(size = fontsizes[1])
                ),
                tickfont = list(size = fontsizes[2])
              ),
              zaxis = list(
                title = list(
                  text = (
                    paste(
                      "Component",
                      compZ,
                      "-",
                      explainedvariance[compZ],
                      "% of variance",
                      sep = " "
                    )
                  ),
                  font = list(size = fontsizes[1])
                ),
                tickfont = list(size = fontsizes[2])
              )
            ),
            title = list(text = "PCA Score Plot",
                         font = list(size = fontsizes[3])),
            legend = list(font = list(size = fontsizes[4]))
          )
          
          return(pca_score_plot)
        }
        
      }
      if (treatment == 4) {
        #2 factors to show
        
        if (is.null(compZ)) {
          # Plot 2D, 2 factors.
        #   if(plottype == "Score_ellipse"){
        #     
        #     colcolcol <- paleta_tratamientos[1:length(unique(tratamientocolores))]
        #     if (is.null(compZ)) {# Plot 2D, 1 factor.
        #       
        #       h <- ggplot2::ggplot(as.data.frame(pca_object$x), 
        #                            aes(pca_object$x[, compX],
        #                                pca_object$x[, compY], 
        #                                color = treatmentstable[,3] ))+ geom_point()+ggplot2::theme_minimal()
        #       
        #       i <- h + ggplot2::stat_ellipse(level=confidence, geom = "polygon", alpha=0.3, 
        #                                      aes(fill=treatmentstable[,3]))+ ggplot2::theme(legend.title = element_blank())
        #       j <- i + ggplot2::scale_color_manual(values = colcolcol) + ggplot2::scale_fill_manual(values=colcolcol)
        #       pca_score_plot3 <- plotly::ggplotly(j)
        #       pca_score_plot3 <- plotly::style(pca_score_plot3, showlegend = F, traces = c(5,6,7,8))
        #       pca_score_plot2 <- plotly::layout(
        #         pca_score_plot3,
        #         
        #         xaxis = list(
        #           zeroline=T,
        #           title = list(
        #             text = paste(
        #               "Component",
        #               compX,
        #               "-",
        #               explainedvariance[compX],
        #               "% of variance",
        #               sep = " "
        #             ),
        #             font =list(
        #               size = fontsizes[1]
        #             )
        #           ),
        #           tickfont = list(size =
        #                             fontsizes[2])
        #         ),
        #         yaxis = list(
        #           zeroline=T,
        #           title = list(
        #             text = paste(
        #               "Component",
        #               compY,
        #               "-",
        #               explainedvariance[compY],
        #               "% of variance",
        #               sep = " "
        #             ),
        #             font =
        #               list(size = fontsizes[1])
        #           ),
        #           tickfont = list(size =
        #                             fontsizes[2])
        #         ),
        #         title = list(text = "PCA Score Plot", font =
        #                        list(size = fontsizes[3])),
        #         legend = list(font = list(size =
        #                                     fontsizes[4]))
        #       )
        #       pca_score_plot <- plotly::style(pca_score_plot2, hoverinfo='none', traces = c(5,6,7,8))
        #       
        #       return(pca_score_plot)
        #   } else{
          pca_score_plot <-
            plotly::plot_ly(
              as.data.frame(pca_object$x),
              x = ~ pca_object$x[, compX],
              y = ~ pca_object$x[, compY],
              hoverinfo = "text",
              text = paste("Sample:", rownames(pca_object$x), sep = " "),
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
          pca_score_plot <- plotly::layout(
            pca_score_plot,
            xaxis = list(
              title = list(
                text = (
                  paste(
                    "Component",
                    compX,
                    "-",
                    explainedvariance[compX],
                    "% of variance",
                    sep = " "
                  )
                ),
                font =
                  list(size = fontsizes[1])
              ),
              tickfont = list(size =
                                fontsizes[2])
            ),
            yaxis = list(
              title = list(
                text = (
                  paste(
                    "Component",
                    compY,
                    "-",
                    explainedvariance[compY],
                    "% of variance",
                    sep = " "
                  )
                ),
                font =
                  list(size = fontsizes[1])
              ),
              tickfont = list(size =
                                fontsizes[2])
            ),
            title = list(text = "PCA Score Plot", font =
                           list(size = fontsizes[3])),
            legend = list(font = list(size =
                                        fontsizes[4]))
          )

          return(pca_score_plot)

          }
        }
        
        if (!is.null(compZ)) {
          # Plot 3D, 2 factors.
          if (compZ > ncol(pca_object$x))
            stop(
              paste(
                "\nZ-Axis component do not exist. Only ",
                length(pca_object$x),
                " components were defined in PCA analysis",
                sep = ""
              )
            )
          
          pca_score_plot <-
            plotly::plot_ly(
              as.data.frame(pca_object$x),
              x = ~ pca_object$x[, compX],
              y = ~ pca_object$x[, compY],
              z = ~ pca_object$x[, compZ],
              hoverinfo = "text",
              text = paste("Sample:", rownames(pca_object$x), sep =
                             " "),
              showlegend = TRUE,
              type = "scatter3d",
              mode = "markers",
              color = factor(treatmentstable[, 1], levels = unique(treatmentstable[, 1])),
              colors = paleta_tratamientos[1:length(unique(treatmentstable[, 1]))],
              symbol = factor(treatmentstable[, 2], levels =
                                unique(treatmentstable[, 2])),
              symbols = simbolos[1:length(unique(treatmentstable[, 2]))],
              marker = list(size = 11)
            )
          pca_score_plot <- plotly::layout(
            pca_score_plot,
            scene = list(
              xaxis = list(
                title = list(
                  text = (
                    paste(
                      "Component",
                      compX,
                      "-",
                      explainedvariance[compX],
                      "% of variance",
                      sep = " "
                    )
                  ),
                  font = list(size = fontsizes[1])
                ),
                tickfont = list(size = fontsizes[2])
              ),
              yaxis = list(
                title = list(
                  text = (
                    paste(
                      "Component",
                      compY,
                      "-",
                      explainedvariance[compY],
                      "% of variance",
                      sep = " "
                    )
                  ),
                  font = list(size = fontsizes[1])
                ),
                tickfont = list(size = fontsizes[2])
              ),
              zaxis = list(
                title = list(
                  text = (
                    paste(
                      "Component",
                      compZ,
                      "-",
                      explainedvariance[compZ],
                      "% of variance",
                      sep = " "
                    )
                  ),
                  font = list(size = fontsizes[1])
                ),
                tickfont = list(size = fontsizes[2])
              )
            ),
            title = list(text = "PCA Score Plot", font = list(size = fontsizes[3])),
            legend = list(font = list(size = fontsizes[4]))
          )
          
          return(pca_score_plot)
          
        }
      }
   
    # Scree plot ----
    if (plottype == "Scree") {
      explainedvariance <-
        pca_object$sdev ^ 2  #Explained variance for each treatment
      screeplotdata <-
        data.frame(
          PCs = paste("PC", c(1:length(
            explainedvariance
          )), sep = ""),
          ExplVar = round(explainedvariance * 100 / sum(explainedvariance), 2),
          Cumulat = round(cumsum(explainedvariance) * 100 / sum(explainedvariance), 2),
          stringsAsFactors = FALSE
        ) #Creamos tabla de datos con % varianza y % var acumulado
      screeplotdata$PCs <-
        factor(screeplotdata$PCs, levels = unique(screeplotdata$PCs)[order(screeplotdata$ExplVar, decreasing = TRUE)]) # Ordenamos los PCs segun valor
      
      pca_screeplot <- plotly::plot_ly(screeplotdata)
      pca_screeplot <-
        plotly::add_trace(
          pca_screeplot,
          x =  ~ PCs,
          y =  ~ ExplVar,
          type = "bar",
          name = "% Explained Variance",
          marker = list(color = '#C9EFF9')
        )
      pca_screeplot <-
        plotly::add_trace(
          pca_screeplot,
          x =  ~ PCs,
          y =  ~ Cumulat,
          type = "scatter",
          mode = "lines+markers",
          yaxis = "y2",
          name = "Cumulative Proportion",
          line = list(color = '#45171D')
        )
      pca_screeplot <- plotly::layout(
        pca_screeplot,
        title = list(text = "PCA Scree Plot", font =
                       list(size = fontsizes[3])),
        xaxis = list(
          title = "",
          titlefont = list(size = fontsizes[1]),
          tickfont = list(size = fontsizes[2])
        ),
        yaxis = list(
          side = 'left',
          title = "% of Explained Variance",
          showgrid = F,
          zeroline = TRUE,
          titlefont = list(size = fontsizes[1]),
          tickfont = list(size = fontsizes[2])
        ),
        yaxis2 = list(
          automargin = TRUE,
          side = 'right',
          overlaying = "y",
          title = "Cumulative",
          showgrid = FALSE,
          zeroline = FALSE,
          titlefont = list(size = fontsizes[1]),
          tickfont = list(size = fontsizes[2])
        ),
        legend = list(x = 0.7, y = 0.1),
        font = list(size = fontsizes[4])
      )
      
      return(pca_screeplot)
      
    }
    
    # Var plot ----
    if (plottype == "Var") {
    
      if (useannot == FALSE) {
        pca_VarPlot <- plotly::plot_ly(
          as.data.frame(pca_object$rotation),
          x = ~ pca_object$rotation[, compX],
          y = ~ pca_object$rotation[, compY],
          text = paste("ID:", rownames(pca_object$rotation)),
          type = "scatter",
          mode = "markers",
          marker = list(size = 6, opacity = 1),
          colors = paleta_niveles[1:length(datasetnames)]
        )
        
        pca_VarPlot <- plotly::layout(
          pca_VarPlot,
          xaxis = list(
            title = list(
              text = (
                paste(
                  "Component",
                  compX,
                  "-",
                  explainedvariance[compX],
                  "% of variance",
                  sep = " "
                )
              ),
              font = list(size = fontsizes[1])
              ),
            tickfont = list(size = fontsizes[2])
          ),
          yaxis = list(
            title = list(
              text = (
                paste(
                  "Component",
                  compY,
                  "-",
                  explainedvariance[compY],
                  "% of variance",
                  sep = " "
                )
              ),
              font = list(size = fontsizes[1])
              ),
            tickfont = list(size = fontsizes[2])
          ),
          title = list(text = "PCA Variable Plot", 
                       font = list(
                         size = fontsizes[3])),
          legend = list(font = list(size = fontsizes[4]))
        )
        
        
        return(pca_VarPlot)
      }
      if (useannot == TRUE) {
        #Annotations
        if (!is.null(compZ))
          cat("\nVarplot can only plot 2 components. Only plotting X and Y.")
        
        annotations2 <- as.data.frame(row.names(pca_object$rotation))
        colnames(annotations2) <- "IDENTIFIER"
        annotations <-
          dplyr::left_join(annotations2, as.data.frame(annotations), by = "IDENTIFIER")
        annotations <-
          as.matrix(annotations) #required to maintain order along with rownames, colors and descriptors
        rm(annotations2)
      
        #chapuza para las anotaciones. eliminamos duplicados pa porsi. en base al ID
        annotations <- annotations[!duplicated(annotations[, 1]), ]
        ####FLAG!!!!!!!!!!!!!!!
        #Sacamos los niveles del mapman, y las etiquetas
        mapmanlevels <-
          as.numeric(unique(annotations[order(as.numeric(annotations[, 3])), 3]))
        mapmannames <-
          annotations[order(as.numeric(annotations[, 3])), ]
        mapmannames <- as.vector(unique(mapmannames[, 4]))
        
        pca_VarPlot <-
          plotly::plot_ly(
            as.data.frame(pca_object$rotation),
            x = ~ pca_object$rotation[, compX],
            y = ~ pca_object$rotation[, compY],
            text = paste(
              "ID:",
              rownames(pca_object$rotation),
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
        
        pca_VarPlot <- plotly::layout(
          pca_VarPlot,
          
          xaxis = list(
            title = list(
              text = paste(
                "Component",
                compX,
                "-",
                explainedvariance[compX],
                "% of variance",
                sep = " "
              ),
              font = list(size = fontsizes[1])
            ),
            tickfont = list(size = fontsizes[2])
          ),
          
          yaxis = list(
            title = list(
              text = paste(
                "Component",
                compY,
                "-",
                explainedvariance[compY],
                "% of variance",
                sep = " "
              ),
              font = list(size = fontsizes[1])
            ),
            tickfont = list(size = fontsizes[2])
          ),
          
          title = list(
            text = "PCA Variable Plot", 
            font = list(size = fontsizes[3])
            ),
          
          legend = list(
            font = list(size = fontsizes[4]))
        )
        
        return(pca_VarPlot)
      }
    }
    
    # Biplot ----
    if (plottype == "Biplot") {
        # Annotations ----
      if (useannot == TRUE) {
        annotations2 <- as.data.frame(row.names(pca_object$rotation))
        colnames(annotations2) <- "IDENTIFIER"
        annotations <-
          dplyr::left_join(annotations2, as.data.frame(annotations), by = "IDENTIFIER")
        annotations <- as.matrix(annotations)
        rm(annotations2)
      }
      
        # Biplot normalizations ----
      scores <- pca_object$x  # Score matrix
      scale <-
        0.9 #scaling factor to merge scores and rotations in one figure default 0.9 /recommended range: 0.9 - 0.7
      lam <- pca_object$sdev[c(compX, compY)] # sqrt(eigenvalues)
      n <- nrow(pca_object$x)
      lam <- lam * sqrt(n) # lam scaling, see below
      # Los scores de prcomp tienen varianza igual al autovalor. Si los dividimos por la raiz cuadrada del autovalor escalaremos por unidad de varianza.
      # Si queremos escalar por la suma de cuadrados, hay que multiplicarlo por la raiz de n
      
      x <- scores[, c(compX, compY)] / lam  # score scaling
      y <- pca_object$rotation[, c(compX, compY)] * lam # rotation scaling
      
      unsigned.range <-
        function(x)
          c(-abs(min(x, na.rm = TRUE)), abs(max(x, na.rm = TRUE))) # sacamos cuales son los rangos para cada eje de scores y loadings
      rangx1 <- unsigned.range(x[, 1L])
      rangx2 <- unsigned.range(x[, 2L])
      rangy1 <- unsigned.range(y[, 1L])
      rangy2 <- unsigned.range(y[, 2L])
      
      ratioejex <- rangx1 / rangy1
      ratioejey <- rangx2 / rangy2
      
      ratio <- max(ratioejex, ratioejey)
      y = y * ratio * scale
      
        # Biplot figure ----
      if (treatment != 4) {
        if (treatment == 1)
          tratamientocolores <- treatmentstable[, 1]
        if (treatment == 2)
          tratamientocolores <- treatmentstable[, 2]
        if (treatment == 3)
          tratamientocolores <- treatmentstable[, 3]
        
        tratamientocolores <-
          factor(tratamientocolores, levels = unique(tratamientocolores))
        
        if (useannot == FALSE) {
          #Mapman levels
          paletaintegrada <-
            c("#466fc7", paleta_tratamientos[1:length(unique(tratamientocolores))])
          
          pca_biplot <- plotly::plot_ly()
          pca_biplot <- plotly::add_markers(
            pca_biplot,
            x = ~ y[, 1],
            y = ~ y[, 2],
            text = paste("ID:", rownames(pca_object$rotation)),
            type = "scatter",
            mode = "markers",
            marker = list(size = 6, opacity = 1),
            color = as.factor("Loadings"),
            colors = paletaintegrada
          )
          pca_biplot <- plotly::add_markers(
            pca_biplot,
            x = ~ x[, 1],
            y = ~ x[, 2],
            inherit = FALSE,
            hoverinfo = "text",
            text = paste("Sample:", rownames(x), sep =
                           " "),
            showlegend = TRUE,
            type = "scatter",
            mode = "markers",
            color = tratamientocolores,
            colors = paletaintegrada[2:length(paletaintegrada)],
            marker = list(size = 11),
            yaxis = "y",
            xaxis = "x"
          )
          
          pca_biplot <- plotly::layout(
            pca_biplot,
            xaxis = list(
              title = list(
                text = paste(
                  "Component",
                  compX,
                  "-",
                  explainedvariance[compX],
                  "% of variance",
                  sep = " "
                ),
                font = list(
                  size = fontsizes[1])
              ),
              tickfont = list(size = fontsizes[2])
            ),
            
            yaxis = list(
              title =  list(
                text = paste(
                  "Component",
                  compY,
                  "-",
                  explainedvariance[compY],
                  "% of variance",
                  sep = " "
                ),
                font = list(
                  size = fontsizes[1])
              ),
              tickfont = list(size = fontsizes[2])
            ),
            title = list(text = "PCA Biplot", font =
                           list(size = fontsizes[3])),
            legend = list(font = list(size = fontsizes[4]))
          )
          
          return(pca_biplot)
        }
        if (useannot == TRUE) {
          #Annotations fixes ... remove suplicates
          annotations <- annotations[!duplicated(annotations[, 1]), ]
          
          # Mapman levels merge and taggs
          mapmanlevels <-
            as.numeric(unique(annotations[order(as.numeric(annotations[, 3])), 3]))
          mapmannames <-
            annotations[order(as.numeric(annotations[, 3])), ]
          mapmannames <- as.vector(unique(mapmannames[, 4]))
          paletaintegrada <-
            c(paleta_mapman[mapmanlevels], paleta_tratamientos[1:length(unique(tratamientocolores))])
          
          pca_biplot <- plotly::plot_ly()
          pca_biplot <- plotly::add_markers(
            pca_biplot,
            x = ~ y[, 1],
            y = ~ y[, 2],
            text = paste(
              "ID:",
              rownames(pca_object$rotation),
              "</br> Description: ",
              annotations[, 2],
              "</br> Mercator Bin: ",
              annotations[, 4]
            ),
            type = "scatter",
            mode = "markers",
            marker = list(size = 6, opacity =
                            1),
            color = factor(annotations[, 4], levels =
                             mapmannames),
            colors = paletaintegrada
          )
          pca_biplot <- plotly::add_markers(
            pca_biplot,
            x = ~ x[, 1],
            y = ~ x[, 2],
            inherit = FALSE,
            hoverinfo = "text",
            text = paste("Sample:", rownames(x), sep =
                           " "),
            showlegend = TRUE,
            type = "scatter",
            mode = "markers",
            color = factor(tratamientocolores, levels = c(unique(
              tratamientocolores
            ))),
            colors = paletaintegrada[(length(paletaintegrada) -
                                        length(unique(tratamientocolores))):(length(paletaintegrada))],
            marker = list(size = 11),
            yaxis = "y",
            xaxis = "x"
          )
          
          pca_biplot <- plotly::layout(
            pca_biplot,
            
            xaxis = list(
              title = list(
                text = (paste(
                  "Component",
                  compX,
                  "-",
                  explainedvariance[compX],
                  "% of variance",
                  sep = " "
                )),
                font = list(size = fontsizes[1])
              ),
              tickfont = list(size = fontsizes[2])
            ),
            
            yaxis = list(
              title =  list(
                text = (paste(
                  "Component",
                  compY,
                  "-",
                  explainedvariance[compY],
                  "% of variance",
                  sep = " "
                  )
                ),
                font = list(size = fontsizes[1])
              ),
              tickfont = list(size = fontsizes[2])
            ),
            
            title = list(text = "PCA Biplot", font =
                           list(size = fontsizes[3])),
            legend = list(font = list(size = fontsizes[4]))
          )
          
          return(pca_biplot)
        }
      }
      
      if (treatment == 4) {
        if (useannot == FALSE) {
          paletaintegrada <-
            c("#466fc7", paleta_tratamientos[1:length(unique(treatmentstable[, 1]))])
          
          pca_biplot <- plotly::plot_ly()
          pca_biplot <- plotly::add_markers(
            pca_biplot,
            x = ~ y[, 1],
            y = ~ y[, 2],
            text = paste("ID:", rownames(pca_object$rotation)),
            type = "scatter",
            mode = "markers",
            marker = list(size = 6, opacity =
                            1),
            color = as.factor("Loadings"),
            colors = paletaintegrada
          )
          
          pca_biplot <- plotly::add_trace(
            pca_biplot,
            x = ~ x[, 1],
            y = ~ x[, 2],
            hoverinfo = "text",
            text = paste("Sample:", rownames(x), sep =
                           " "),
            showlegend = TRUE,
            type = "scatter",
            mode = "markers",
            color = factor(treatmentstable[, 1], levels = unique(treatmentstable[, 1])),
            colors = paletaintegrada[2:length(paletaintegrada)],
            symbol = factor(treatmentstable[, 2], levels =
                              unique(treatmentstable[, 2])),
            symbols = simbolos[1:length(unique(treatmentstable[, 2]))],
            marker = list(size = 11),
            inherit = F
          )
          
          pca_biplot <- plotly::layout(
            pca_biplot,
            
            xaxis = list(
              title = list(
                text = (paste(
                  "Component",
                  compX,
                  "-",
                  explainedvariance[compX],
                  "% of variance",
                  sep = " "
                )),
              font = list(
                size = fontsizes[1])
              ),
              tickfont = list(size = fontsizes[2])
            ),
            
            yaxis = list(
              title = list(
                text = (paste(
                  "Component",
                  compY,
                  "-",
                  explainedvariance[compY],
                  "% of variance",
                  sep = " ")
                  ),
                font = list(size = fontsizes[1])
                ),
              tickfont = list(size = fontsizes[2])
            ),
            
            title = list(
              text = "PCA Biplot", 
              font = list(
                size = fontsizes[3])),
            
            legend = list(
              font = list(
                size = fontsizes[4]
                )
              )
          )
          
          return(pca_biplot)
          
        }
        if (useannot == TRUE) {
          annotations <- annotations[!duplicated(annotations[, 1]), ]
          #FLAG
          mapmanlevels <-
            as.numeric(unique(annotations[order(as.numeric(annotations[, 3])), 3]))
          mapmannames <-
            annotations[order(as.numeric(annotations[, 3])), ]
          mapmannames <- as.vector(unique(mapmannames[, 4]))
          paletaintegrada <-
            c(paleta_mapman[mapmanlevels], paleta_tratamientos[1:length(unique(treatmentstable[, 1]))])
          
          pca_biplot <- plotly::plot_ly()
          pca_biplot <- plotly::add_markers(
            pca_biplot,
            x = ~ y[, 1],
            y = ~ y[, 2],
            text = paste(
              "ID:",
              rownames(pca_object$rotation),
              "</br> Description: ",
              annotations[, 2],
              "</br> Mercator Bin: ",
              annotations[, 4]
            ),
            type = "scatter",
            mode = "markers",
            marker = list(size = 6, opacity =
                            1),
            color = factor(annotations[, 4], levels =
                             mapmannames),
            colors = paletaintegrada
          )
          
          pca_biplot <- plotly::add_markers(
            pca_biplot,
            x = ~ x[, 1],
            y = ~ x[, 2],
            inherit = FALSE,
            hoverinfo = "text",
            text = paste("Sample:", rownames(x), sep =" "),
            showlegend = TRUE,
            type = "scatter",
            mode = "markers",
            color = factor(treatmentstable[, 1], levels = c(unique(
              treatmentstable[, 1]
            ))),
            colors = paletaintegrada[(length(paletaintegrada) -
                                        length(unique(treatmentstable[, 1]))):(length(paletaintegrada))],
            marker = list(size = 11),
            yaxis = "y",
            xaxis = "x",
            symbol = factor(treatmentstable[, 2], levels =
                              unique(treatmentstable[, 2])),
            symbols = simbolos[1:length(unique(treatmentstable[, 2]))],
            marker = list(size = 11)
          )
          
          pca_biplot <- plotly::layout(
            pca_biplot,
            xaxis = list(
              title = list(
                text = (paste(
                  "Component",
                  compX,
                  "-",
                  explainedvariance[compX],
                  "% of variance",
                  sep = " "
                  )
                ),
              font = list(size = fontsizes[1])),
              tickfont = list(size = fontsizes[2])
            ),
            
            yaxis = list(
              title = list(
                text = (paste(
                  "Component",
                  compY,
                  "-",
                  explainedvariance[compY],
                  "% of variance",
                  sep = " "
                  )
                ),
              font = list(size = fontsizes[1])
              ),
              tickfont = list(size = fontsizes[2])
            ),
            
            title = list(text = "PCA Biplot", 
                         font = list(size = fontsizes[3])),
            
            legend = list(font = list(size = fontsizes[4]))
          )
          
          
          return(pca_biplot)
        }
      }
      
    }
    # Composite Plot----
    if (plottype == "Composite") {
      if (!is.null(compZ))
        stop("Composite plot can only draw 2D plots. compZ must be null")
      
      p1 <-
        pca_plot(
          pca_list,
          treatment,
          compX,
          compY,
          plottype = "Score",
          useannot = useannot,
          fontsizes = c(11, 10, 16, 12)
        )
      plotly::layout(
        p1,
        annotations = list(
          text = "Score Plot",
          xref = "paper",
          yref = "paper",
          font = list(size = 14),
          margin = 1,
          yanchor = "bottom",
          xanchor = "center",
          align = "center",
          x = 0.5,
          y = 1,
          showarrow = FALSE
        )
      )
      p2 <-
        pca_plot(
          pca_list,
          treatment,
          compX,
          compY,
          plottype = "Var",
          useannot = useannot,
          fontsizes = c(11, 10, 16, 12)
        )
      plotly::layout(
        p2,
        annotations = list(
          text = "Variable Plot",
          xref = "paper",
          yref = "paper",
          font = list(size = 14),
          margin = 1,
          yanchor = "bottom",
          xanchor = "center",
          align = "center",
          x = 0.5,
          y = 1,
          showarrow = FALSE
        )
      )
      p3 <-
        pca_plot(
          pca_list,
          treatment,
          compX,
          compY,
          plottype = "Biplot",
          useannot = useannot,
          fontsizes = c(11, 10, 16, 12)
        )
      plotly::layout(
        p3,
        annotations = list(
          text = "Biplot",
          xref = "paper",
          yref = "paper",
          font = list(size = 14),
          margin = 1,
          yanchor = "bottom",
          xanchor = "center",
          align = "center",
          x = 0.5,
          y = 1,
          showarrow = FALSE
        )
      )
      #p4 <- pca_plot(pca_list,treatment, compX, compY, plottype = "Scree")layout(annotations=list(text = "Scree Plot",xref = "paper",yref = "paper",
      #yanchor = "bottom",xanchor = "center",align = "center",x = 0.5,y = 1,showarrow = FALSE))
      #mcoa_ScorePlot <- subplot(subplot(p1,p2),subplot(p3,p4),nrows=2,shareX = F,shareY = F,margin = 0.04)
      pca_CompositeplotPlot <-
        plotly::subplot(
          plotly::subplot(p1, p2, titleX = TRUE, titleY = TRUE),
          plotly::subplot(p3, titleX = TRUE, titleY = TRUE),
          nrows = 2,
          shareX = F,
          shareY = F,
          titleX = TRUE,
          titleY = TRUE
        )
      pca_CompositeplotPlot <- plotly::layout(
        pca_CompositeplotPlot,
        title = list(text = "PCA Analysis", font = list(size =
                                                          fontsizes[3])),
        autosize = TRUE
      )
      return(pca_CompositeplotPlot)
    }
    # All ----
    if (plottype == "All") {
      #Recursive calls and export_plot() that's it
      if (!is.null(compZ))
        stop("Please select only two dimensions for plotting all the figures")
      cat("EXPORTING PCA ANALYSIS PLOTS\nA single core will be used. It may take a while...\n ")
      cat("Generating and saving Score plot. Step 1/5\n")
      export_plot(
        pca_plot(
          pca_list,
          treatment = treatment,
          compX = compX,
          compY = compY,
          plottype = "Score",
          useannot = useannot
        ),
        "pca_score.pdf"
      )
      cat("Generating and saving Scree plot. Step 2/5\n")
      export_plot(
        pca_plot(
          pca_list,
          treatment = treatment,
          compX = compX,
          compY = compY,
          plottype = "Scree",
          useannot = useannot
        ),
        "pca_scree.pdf"
      )
      cat("Generating and saving Var plot. Step 3/5\n")
      export_plot(
        pca_plot(
          pca_list,
          treatment = treatment,
          compX = compX,
          compY = compY,
          plottype = "Var",
          useannot = useannot
        ),
        "pca_var.pdf"
      )
      cat("Generating and saving Biplot. Step 4/5\n")
      export_plot(
        pca_plot(
          pca_list,
          treatment = treatment,
          compX = compX,
          compY = compY,
          plottype = "Biplot",
          useannot = useannot
        ),
        "pca_biplot.pdf"
      )
      cat("Generating and saving Composite plot. Step 5/5\n")
      export_plot(
        pca_plot(
          pca_list,
          treatment = treatment,
          compX = compX,
          compY = compY,
          plottype = "Composite",
          useannot = useannot
        ),
        "pca_composite.pdf"
      )
    }
  }
