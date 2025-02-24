

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

#' @name ica_plot
#' @title A function to plot Independent Components Analysis (ICA)
#' @description This function makes different plots of ica analysis
#' @usage ica_plot (ica_list, treatment = 1, compX = 1, compY = 2,
#' compZ = NULL, plottype = "Composite", useannot = FALSE,
#' fontsizes = c(14,10,16,12), fortopscoring = c(1,20,"abs"),
#' confidence = 0.9)
#' @param ica_list ica object to plot, class "icaanalysis"
#' @param treatment numeric value to determine the colors and symbols to be displayed in the plot, used in Score plot and Topscring plot. Available optiones are:
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
#' \item Topscoring, top scoring variables in the selected (IC) independent component
#' \item Var, variable plot, each dot represent a variable.
#' \item Composite: combination of the two most representative plots (Score and Var)
#' \item All all aforementioned plots will be exported as .pdf files in your working directory.
#' }
#' @param useannot logical, indicating if annotation should be used or not, please note the previous ica analysis should include it
#' @param fortopscoring vector indicating the number of the component to be analyzed, the top n loadings to plot and the mode to compute "abs" for absolute value and "posneg" for positive and negative loading analyzed separately
#' @param fontsizes vector containing the font sizes to use in the plots
#' @param confidence numeric value between 0 and 1 for confidence ellipses in Score_Ellipse plot
#' @return ica plot
#' @author Luis Valledor and Laura Lamelas

#' @importFrom plotly plot_ly layout add_trace subplot style ggplotly
#' @importFrom dplyr left_join
#' @importFrom ggplot2 ggplot theme scale_color_manual scale_fill_manual stat_ellipse geom_point theme_minimal aes
#' @importFrom stats na.omit

#' @export

ica_plot <- function(ica_list,
                     treatment = 1,
                     compX = 1,
                     compY = 2,
                     compZ = NULL,
                     plottype = "Composite",
                     useannot = FALSE,
                     fontsizes = c(14, 10, 16, 12),
                     fortopscoring = c(1, 20, "abs"),
                     confidence = 0.9) {

  # Flags ----
  explainedvariance <- NULL  
  # Checks ----
  if (!inherits(ica_list,"icaanalysis"))
    stop("an ICA object is expected, please run ica_analysis first")
  if (!is.list(ica_list[[1]]))
    stop("ICA analysis object seems to be wrong. Please re-run ica_analysis")
  if (plottype %in% c("Composite", "Score","Score_Ellipse", "Scree", "Var", "Topscoring", "All") == FALSE)
    stop("Please select a valid plot type")
  if (useannot == T && is.null(ica_list$annotations)) {
    cat(
      paste(
        "\nAnnotation matrix not provided. Considering useannot=F\n",
        "Please check/repeat ica_analysis() including an annotation matrix\n",
        "Read help or follow the tutorial for more information.\n",
        sep = ""
      )
    )
    useannot = FALSE
  }
  
  # Required objects ----
  ica_object <- ica_list$ica
  treatmentstable <- ica_list$treatments
  annotations <- ica_list$annotations
  datasetnames <- ica_list$datasetnames
  treatmentnames <- colnames(treatmentstable)
  originaldata <- ica_list$originaldata
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
          ica_list,
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
          ica_list,
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
  # Score plot ----
  if (plottype %in% c("Score", "Score_Ellipse")) {
    if (compX > length(ica_object$vafs))
      stop(
        paste(
          "\nX-Axis component do not exist. Only ",
          length(ica_object$vafs),
          " components were defined in ICA analysis",
          sep = ""
        )
      )
    if (compY > length(ica_object$vafs))
      stop(
        paste(
          "\nY-Axis component do not exist. Only ",
          length(ica_object$vafs),
          " components were defined in ICA analysis",
          sep = ""
        )
      )
    
    
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
      
      
      if(plottype == "Score_Ellipse"){
        replicates <- nrow(treatmentstable)/length(unique(tratamientocolores))
        if (replicates < 4){plottype = "Score"}
        if (!is.null(compZ)){plottype = "Score"}
        
      }
      if(plottype == "Score_Ellipse"){
        colcolcol <- paleta_tratamientos[1:length(unique(tratamientocolores))]
        if (is.null(compZ)) {# Plot 2D, 1 factor.
          
          h <- ggplot2::ggplot(as.data.frame(ica_object$S), 
                               ggplot2::aes(ica_object$S[, compX],
                                            ica_object$S[, compY], 
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
          
          ica_score_plot3 <- plotly::ggplotly(j)
          
          ica_score_plot2 <- plotly::layout(
            ica_score_plot3,
            
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
            title = list(text = "ICA Score Plot", font =
                           list(size = fontsizes[3])),
            legend = list(font = list(size =
                                        fontsizes[4]))
          )
          
          #retocamos trazas elipses
          ica_score_plot2_t <- ica_score_plot2
          ellips_text <- paste(confidence*100,"% confidence","<br>interval; ",unique(tratamientocolores))
          
          for(i in (1:length(unique(colcolcol)))){
            ica_score_plot2_t <- plotly::style(ica_score_plot2_t, 
                                               hoverinfo = "text", text = ellips_text[i], traces = i,name = ellips_text[i], showlegend=T)
          }
          #retocamos trazas samples
          ica_score_plot_tp <- ica_score_plot2_t
          track <- c((1+length(unique(colcolcol))):(length(unique(colcolcol))*2)) 
          ptra <- as.data.frame(cbind(rownames(ica_object$S),treatmentstable[,treatment]))
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
            ica_score_plot_tp <- plotly::style(ica_score_plot_tp, 
                                               hoverinfo = "text",
                                               text = paste("Sample: ",stats::na.omit(b[,i]), sep=""),
                                               showlegend = TRUE, 
                                               name = colnames(b)[i],
                                               traces = track[i] )}
          return(ica_score_plot_tp)}
      }
      
      
      else {# Plot 2D, 1 factor.
        ica_score_plot <-
          plotly::plot_ly(
            as.data.frame(ica_object$S),
            x = ~ ica_object$S[, compX],
            y = ~ ica_object$S[, compY],
            hoverinfo = "text",
            text = paste("Sample:", rownames(ica_object$S), sep =" "),
            showlegend = TRUE,
            type = "scatter",
            mode = "markers",
            color = tratamientocolores,
            colors = paleta_tratamientos[1:length(unique(tratamientocolores))],
            #symbol = factor(rep(c(1:nrow(mcoin_object$mcoa$lambda)),each=nrow(mcoin_object$mcoa$SynVar)), levels=datasetnames),
            symbols = simbolos[1:length(unique(tratamientocolores))],
            marker = list(size = 11)
          )
        ica_score_plot <- plotly::layout(
            ica_score_plot,
            xaxis = list(
              title = list(
                text = paste(
                "Component",
                compX,
                "-",
                round(ica_object$vafs * 100, 2)[compX],
                "% of variance",
                sep = " "
              ),
              font = 
                list(size = fontsizes[1])
              ),
              tickfont = list(size = fontsizes[2])
            ),
            yaxis = list(
              title = list (text = paste(
                "Component",
                compY,
                "-",
                round(ica_object$vafs * 100, 2)[compY],
                "% of variance",
                sep = " "
              ),
                font = list(size = fontsizes[1])),
              tickfont = list(size = fontsizes[2])
            ),
            title = list(text="ICA ScorePlot",font = list(size = fontsizes[3])),
            legend = list(font = list(size =
                                        fontsizes[4]))
          )
        
        return(ica_score_plot)
      }
      
      if (!is.null(compZ)) {
        # Plot 3D, 1 factor.
        if (compZ > length(ica_object$vafs))
          stop(
            paste(
              "\nZ-Axis component do not exist. Only ",
              length(ica_object$vafs),
              " components were defined in ICA analysis",
              sep = ""
            )
          )
        ica_score_plot <-
          plotly::plot_ly(
            as.data.frame(ica_object$S),
            x = ~ ica_object$S[, compX],
            y = ~ ica_object$S[, compY],
            z = ~ ica_object$S[, compZ],
            hoverinfo = "text",
            text = paste("Sample:", rownames(ica_object$S), sep =
                           " "),
            showlegend = TRUE,
            type = "scatter3d",
            mode = "markers",
            color = tratamientocolores,
            colors = paleta_tratamientos[1:length(unique(tratamientocolores))],
            #symbol = tratamientocolores
            #symbols = simbolos[1:length(unique(tratamientocolores))],
            marker = list(size = 11)
          )
        
        ica_score_plot <-
          plotly::layout(
            p = ica_score_plot,
            scene = list(
              xaxis = list(
                title = list( text = paste(
                  "Component",
                  compX,
                  "-",
                  round(ica_object$vafs * 100, 2)[compX],
                  "% of variance",
                  sep = " "
                ),
                font = list(size = fontsizes[1])),
                tickfont = list(size = fontsizes[2])
              ),
              yaxis = list(
                title = list(
                  text = paste(
                    "Component",
                    compY,
                    "-",
                    round(ica_object$vafs * 100, 2)[compY],
                    "% of variance",
                    sep = " "
                  ),
                  font = list(size = fontsizes[1])
                ),
                tickfont = list(size = fontsizes[2])
              ),
              zaxis = list(
                title = list(
                  text = paste(
                    "Component",
                    compZ,
                    "-",
                    round(ica_object$vafs * 100, 2)[compZ],
                    "% of variance",
                    sep = " "
                  ),
                  font = list(size = fontsizes[1])
                ),
                tickfont = list(size = fontsizes[2])
              )
            ),
            title = list (text = "ICA ScorePlot",
            font = list(size = fontsizes[3])),
            legend = list(font = list(size =
                                        fontsizes[4]))
          )
        
        return(ica_score_plot)
      }
      
    }
    if (treatment == 4) {
      # 2 factors to show
      
      if (is.null(compZ)) {
        # Plot 2D, 2 factores.
        ica_score_plot <-
          plotly::plot_ly(
            as.data.frame(ica_object$S),
            x = ~ ica_object$S[, compX],
            y = ~ ica_object$S[, compY],
            hoverinfo = "text",
            text = paste("Sample:", rownames(ica_object$S), sep =
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
        ica_score_plot <- plotly::layout(
            ica_score_plot,
            xaxis = list(
              title = list( 
                text = paste(
                "Component",
                compX,
                "-",
                round(ica_object$vafs * 100, 2)[compX],
                "% of variance",
                sep = " "
              ),
                font = list(size = fontsizes[1])
              ),
              tickfont = list(size = fontsizes[2])
            ),
            yaxis = list(
              title = list(
                text = 
                  paste(
                    "Component",
                    compY,
                    "-",
                    round(ica_object$vafs * 100, 2)[compY],
                    "% of variance",
                    sep = " "
                  ),
                  font = list(size = fontsizes[1])
              ),
              tickfont = list(size = fontsizes[2])
            ),
            title = list(text ="ICA ScorePlot",
            font = list(size = fontsizes[3])),
            legend = list(font = list(size =
                                        fontsizes[4]))
        )
        return(ica_score_plot)
      }
      
      if (!is.null(compZ)) {
        # Plot 3D, 2 factores.
        if (compZ > length(ica_object$vafs))
          stop(
            paste(
              "\nZ-Axis component do not exist. Only ",
              length(ica_object$vafs),
              " components were defined in ICA analysis",
              sep = ""
            )
          )
        
        ica_score_plot <-
          plotly::plot_ly(
            as.data.frame(ica_object$S),
            x = ~ ica_object$S[, compX],
            y = ~ ica_object$S[, compY],
            z = ~ ica_object$S[, compZ],
            hoverinfo = "text",
            text = paste("Sample:", rownames(ica_object$S), sep =
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
        ica_score_plot <-
          plotly::layout(
            p = ica_score_plot,
            scene = list(
              xaxis = list(
                title = list(
                  text=(
                    paste(
                      "Component",
                      compX,
                      "-",
                      round(ica_object$vafs * 100, 2)[compX],
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
                      round(ica_object$vafs * 100, 2)[compY],
                      "% of variance",
                      sep = " "
                    )
                  ),
                font = list(size = fontsizes[1])),
                tickfont = list(size = fontsizes[2])
              ),
              zaxis = list(
                title = list(
                  text = (
                    paste(
                      "Component",
                      compZ,
                      "-",
                      round(ica_object$vafs * 100, 2)[compZ],
                      "% of variance",
                      sep = " "
                    )
                  ),
                font = list(size = fontsizes[1])
                ),
                tickfont = list(size = fontsizes[2])
              )
            ),
            title = list(text = "ICA ScorePlot",
                         font = list(size = fontsizes[3])
                         ),
            legend = list(font = list(size =
                                        fontsizes[4]))
          )
        
        return(ica_score_plot)
        
      }
    }
  }
  
  # Scree plot ----
  if (plottype == "Scree") {
    screeplotdata <-
      data.frame(
        ICs = paste("IC", c(1:length(
          ica_object$vafs
        )), sep = ""),
        ExplVar = round(ica_object$vafs * 100, 2),
        Cumulat = round(cumsum(ica_object$vafs * 100), 2),
        stringsAsFactors = FALSE
      ) 
    screeplotdata$ICs <-
      factor(screeplotdata$ICs, levels = unique(screeplotdata$ICs)[order(screeplotdata$ExplVar, decreasing = TRUE)]) 
    
    ica_screeplot <- plotly::plot_ly(screeplotdata)
    ica_screeplot <-
      plotly::add_trace(
        p = ica_screeplot,
        x =  ~ ICs,
        y =  ~ ExplVar,
        type = "bar",
        name = "% Explained Variance",
        marker = list(color = '#C9EFF9')
      )
    ica_screeplot <-
      plotly::add_trace(
        p = ica_screeplot,
        x =  ~ ICs,
        y =  ~ Cumulat,
        type = "scatter",
        mode = "lines+markers",
        yaxis = "y2",
        name = "Cumulative Proportion",
        line = list(color = '#45171D')
      )
    ica_screeplot <-
      plotly::layout(
        ica_screeplot,
        title =list(
          text= "ICA Screeplot",
          font = list(size = fontsizes[3])
        ),
        xaxis = list(
          title = list(font = list(size = fontsizes[1])),
          tickfont = list(size = fontsizes[2])
        ),
        yaxis = list(
          side = 'left',
          title = list(text = "% of Explained Variance",
                       font = list(size = fontsizes[1])),
          showgrid = F,
          zeroline = TRUE,
          tickfont = list(size = fontsizes[2])
        ),
        yaxis2 = list(
          side = 'right',
          overlaying = "y",
          title = list(text = "Cumulative", font = list(size = fontsizes[1])),
          showgrid = FALSE,
          zeroline = FALSE,
          tickfont = list(size = fontsizes[2])
        ),
        legend = list(x = 0.7, y = 0.1),
        font = list(size = fontsizes[4])
      )
    return(ica_screeplot)
    
  }
  # Var plot ----
  if (plottype == "Var") {
    
    if (useannot == FALSE) {
      ica_VarPlot <- plotly::plot_ly(
        as.data.frame(ica_object$W),
        x = ~ ica_object$W[compX,],
        y = ~ ica_object$W[compY,],
        text = paste("ID:", colnames(ica_object$W)),
        type = "scatter",
        mode = "markers",
        marker = list(size = 6, opacity = 1),
        colors = paleta_niveles[1:length(datasetnames)]
      )
      
      ica_VarPlot <- plotly::layout(
        ica_VarPlot,
        xaxis = list(
          title = list(
            text = (
              paste(
                "Component",
                compX,
                "-",
                round(ica_object$vafs * 100, 2)[compX],
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
                round(ica_object$vafs * 100, 2)[compY],
                "% of variance",
                sep = " "
                )
              ),
            font = list(size = fontsizes[1])
          ),
          tickfont = list(size = fontsizes[2])
        ),
        title = list(
          text ="ICA Estimated Unmixing Signals",
          font = list(
            size = fontsizes[3])
        ),
        legend = list(font = list(size = fontsizes[4]))
      )
      
      
      return(ica_VarPlot)
    }
    if (useannot == TRUE) {
      #Annotations
      if (!is.null(compZ))
        cat("\nVarplot can only plot 2 components. Only plotting X and Y.")
      
      annotations2 <- as.data.frame(colnames(ica_object$W))
      colnames(annotations2) <- "IDENTIFIER"
      annotations <-
        dplyr::left_join(annotations2, as.data.frame(annotations), by = "IDENTIFIER") 
      annotations <-
        as.matrix(annotations) # required!
      rm(annotations2)
      ####FLAG!!!!!!!!!!!!!!!
      #Sacamos los niveles del mapman, y las etiquetas
      mapmanlevels <-
        as.numeric(unique(annotations[order(as.numeric(annotations[, 3])), 3]))
      mapmannames <-
        annotations[order(as.numeric(annotations[, 3])),]
      mapmannames <- as.vector(unique(mapmannames[, 4]))
      
      ica_VarPlot <- plotly::plot_ly(
        as.data.frame(ica_object$W),
        x = ~ ica_object$W[compX,],
        y = ~ ica_object$W[compY,],
        text = paste(
          "ID:",
          colnames(ica_object$W),
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
      
      ica_VarPlot <- plotly::layout(
        ica_VarPlot,
        
        xaxis = list(
          title = list(
            text = paste(
              "Component",
              compX,
              "-",
              round(ica_object$vafs * 100, 2)[compX],
              "% of variance",
              sep = " "
            ),
            font = list(
              size = fontsizes[1])
          ),
          tickfont = list(size = fontsizes[2])
        ),
        
        yaxis = list(
          title = list(
            text = paste(
              "Component",
              compY,
              "-",
              round(ica_object$vafs * 100, 2)[compY],
              "% of variance",
              sep = " "
            ),
            font = list(
              size = fontsizes[1])
          ),
          tickfont = list(size = fontsizes[2])
        ),
        
        title = list(
          text = 
            "ICA Estimated Unmixing Signals",
          font = list(
            size = fontsizes[3])
        ),
        
        legend = list(font = list(size = fontsizes[4]))
      )
      
      return(ica_VarPlot)
    }
  }
  
  # Composite plot ----
  if (plottype == "Composite") {
    if (!is.null(compZ))
      stop("Composite plot can only draw 2D plots. compZ must be null")
    
    p1 <-
      ica_plot(ica_list, treatment, compX, compY, plottype = "Score")
    p2 <-
      ica_plot(ica_list,
               treatment,
               compX,
               compY,
               plottype = "Var",
               useannot = useannot)
    ica_composite <-
      plotly::subplot(
        p1,
        p2,
        nrows = 2,
        shareX = F,
        shareY = F,
        margin = 0.04
      )
    ica_composite <-
      plotly::layout(ica_composite, title = "ICA Analysis")
    return(ica_composite)
  }
  # All ----
  if (plottype == "All") {
    #Aqui habra que crear una llamada recursiva a todos los graficos posibles y devolver una lista.
    #tambien recursivo pa los distintos niveles.
    if (!is.null(compZ))
      stop("Please select only two dimensions for plotting all figures")
    cat("EXPORTING ICA ANALYSIS PLOTS\nA single core will be used. It may take a while...\n ")
    cat("Generating and saving Score plot. Step 1/4\n")
    export_plot(
      ica_plot(
        ica_list,
        treatment = treatment,
        compX = compX,
        compY = compY,
        plottype = "Score",
        useannot = useannot
      ),
      "ica_score.pdf"
    )
    cat("Generating and saving Scree plot. Step 2/4\n")
    export_plot(
      ica_plot(
        ica_list,
        treatment = treatment,
        compX = compX,
        compY = compY,
        plottype = "Scree",
        useannot = useannot
      ),
      "ica_scree.pdf"
    )
    cat("Generating and saving Var plot. Step 3/4\n")
    export_plot(
      ica_plot(
        ica_list,
        treatment = treatment,
        compX = compX,
        compY = compY,
        plottype = "Var",
        useannot = useannot
      ),
      "ica_var.pdf"
    )
    cat("Generating and saving Composite plot. Step 4/4\n")
    export_plot(
      ica_plot(
        ica_list,
        treatment = treatment,
        compX = compX,
        compY = compY,
        plottype = "Composite",
        useannot = useannot
      ),
      "ica_composite.pdf"
    )
  }
}
