

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


#' @name bda_plot
#' @title A function to plot DIABLO based analysis
#' @description This function depicts the results obtained from a DIABLO analysis
#' @usage bda_plot(bda_list, treatment = 1, plottype = "Distance",
#' variablespace = "xy", useannot = FALSE, fontsizes = c(14,10,16,12),
#' cutoff = 0.5, fortopscoring = c(1,20,"abs"), confidence = 0.9)
#' @param bda_list diabloanalysis class object containing the diablo analysis of the data
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
#' \item Distance_Ellipse; distance plot with confidence ellipses by tretment
#' \item Topscoring; most representantive loadings according to Distance plot
#' \item Var; variable plot, each dot represent a variable. select dataset and annotation
#' \item Cim;
#' \item Network; network
#' \item Circos; circos plot
#' \item All; all of above
#' }
#' @param variablespace c("x","y","xy")
#' @param useannot logical, indicating if annotation should be used or not, for annotations use the splsda analysis should include it
#' @param fontsizes font sizes for plot elements
#' @param cutoff cutoff value for network prunning
#' @param confidence numeric value between 0 and 1 for confidence ellipses in Distance_Ellipse plot
#' @param fortopscoring vector indicating the number of the component to be analyzed, the top n loadings to plot and the mode to compute "abs" for absolute value and "posneg" for positive and negative loading analyzed separately
#' @return plotly type plot
#' @author Luis Valledor and Laura Lamelas
#' @export
#' @importFrom grDevices dev.cur dev.off
#' @importFrom plotly subplot layout plot_ly 
#' @importFrom mixOmics plotVar cimDiablo network circosPlot
#' @importFrom dplyr left_join


bda_plot <-
  function(bda_list, treatment = 1, plottype = "Distance", variablespace = "xy", useannot = FALSE, 
           fontsizes = c(14, 10, 16, 12), cutoff = 0.5, fortopscoring = c(1, 20, "abs"), confidence = 0.9) {
    # Checks ----
   
    if(!inherits(bda_list, "bdaanalysis")) stop("bda analysis object is expected, please run bda_analysis() first")
    if (plottype %in% c("Distance", "Distance_Ellipse", "Network", "All", "Topscoring", 
                         "Cim", "Circos") == FALSE) 
      stop("Please select a valid plot type")
    if (variablespace %in% c("x", "y", "xy") == FALSE) 
      stop("Select an adequate variable representation space. x, y, or xy ")
    if (useannot == T && is.null(bda_list$annotations)) {
      texts_wizard("\nAnnotation matrix not provided. Considering useannot=F\n")
      texts_wizard("Please check/repeat bda_analysis() including an annotation matrix\n")
      texts_wizard("Read help or follow the tutorial for more information.\n")
      useannot = FALSE
    }
    #Reset graphical device  
    if (grDevices::dev.cur() > 1) {
      grDevices::dev.off()
    }
    #options(warn = -1)
    # Initial objects ----
    compZ <- NULL
    bda_object <- bda_list$bda
    treatmentstable <- bda_list$treatments
    annotations <- bda_list$annotations
    datasetnames <- bda_list$datasetnames
    treatmentnames <- colnames(treatmentstable)
    originaldata <- bda_list$originaldata
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
    
    
    # Top scoring plot ----
    if (plottype == "Topscoring") {
      if (useannot == T) {
        topsco_plot <-
          topscoring(
            bda_list,
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
            bda_list,
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
      
      a <- list()
      for (i in 1:length(bda_object$variates)) {
        a[[i]] <-
          distanceplotdiablo(
            bda_object$variates[[i]],
            bda_object$explained_variance[[i]],
            treatment,
            treatmentstable,
            names(bda_object$explained_variance)[i], confidence, plottype = plottype
          )
      }
      diablo_distance_plot <- plotly::subplot(a, nrows = 2, titleX = T, titleY = T, margin = 0.05)
      diablo_distance_plot <- plotly::layout(
        diablo_distance_plot,
        title = list(
          text = paste("Distance Plots:", paste(
            names(bda_object$explained_variance), collapse = ", "
            ), sep = " "),
          font = list(size = fontsizes[3])
        ),
        legend = list(font = list(size = fontsizes[4]))
      )
      return(diablo_distance_plot)
    }
    # Var plot ----
    if (plottype == "Var") {
      vartable <-
        invisible(mixOmics::plotVar(bda_object, plot = F)) #esto es una chapuza. habria que diseñar funcion para sacar esto. de momento uso mixomics
      xaxislabel <- paste("Component 1") #FLAG
      yaxislabel <- paste("Component 2") #FLAG
      
      if (useannot == FALSE) {
        diablo_VarPlot <- plotly::plot_ly(
          as.data.frame(vartable),
          x = ~ vartable[, 1],
          y = ~ vartable[, 2],
          text = paste("ID:", rownames(vartable)),
          type = "scatter",
          mode = "markers",
          marker = list(size = 6, opacity = 1),
          colors = paleta_niveles[1:length(datasetnames)],
          symbol = factor(vartable[, 3], levels =
                            unique(vartable[, 3])),
          symbols = simbolos[1:length(unique(vartable[, 3]))]
        )
        
        diablo_VarPlot <-
          plotly::layout(
            diablo_VarPlot,
            xaxis = list(
              title = list(text = xaxislabel,
                           font = list(size = fontsizes[1])),
              tickfont = list(size = fontsizes[2])
            ),
            yaxis = list(
              title = list(text = yaxislabel,
                           font = list(size = fontsizes[1])),
              tickfont = list(size = fontsizes[2])
            ),
            title = list(text = "Diablo Variable Plot",
                         font = list(size = fontsizes[3])),
            legend = list(font = list(size = fontsizes[4]))
          )
        
        
        return(diablo_VarPlot)
      }
      if (useannot == TRUE) {
        #Generamos anotaciones
        
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
        
        diablo_VarPlot <- plotly::plot_ly(
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
          symbol = factor(vartable[, 3], levels =
                            unique(vartable[, 3])),
          symbols = simbolos[1:length(unique(vartable[, 3]))]
        )
        
        diablo_VarPlot <-
          plotly::layout(
            diablo_VarPlot,
            xaxis = list(
              title = list(text = xaxislabel,
                           font = list(size = fontsizes[1])),
              tickfont = list(size = fontsizes[2])
            ),
            yaxis = list(
              title = list(text = yaxislabel,
                           font = list(size = fontsizes[1])),
              tickfont = list(size = fontsizes[2])
            ),
            title = list(text = "Diablo Variable Plot",
                         font = list(size = fontsizes[3])),
            legend = list(font = list(size = fontsizes[4]))
          )
        
        return(diablo_VarPlot)
      }
    }
    
    # Cim plot ----
    if (plottype == "Cim") {
      col.index.tr <- unique(as.numeric(factor(treatmentstable[,treatment])))
      tr.cols <- paleta_tratamientos[1:length(unique(treatmentstable[,treatment]))]
      #aux.lvls <- bda_object$names$blocks[1:(length(bda_object$names$blocks)-1)]
      #col.index.lvl <- unique(as.numeric(factor(aux.lvls)))
      #lvl.cols <- paleta_niveles[1:length(col.index.lvl)]
      if (useannot == T) {
        bda_object <- bda_list[[5]]
        suppressMessages(mixOmics::cimDiablo(
          bda_object,
          color.blocks = paleta_niveles[1:(length(bda_object$names$blocks)-1)],
          color.Y = tr.cols[col.index.tr],
          legend.position = 0,
          size.legend = 1,
          margins = c(8, 8)
        ))
        nam <- names(bda_object)
        bda_object[[length(bda_object)+1]] <- tr.cols[col.index.tr]
        bda_object[[length(bda_object)+1]] <- paleta_niveles[1:(length(bda_object$names$blocks)-1)]
        names(bda_object) <- c(nam, "tr.cols","lvl.cols")
        class(bda_object) <- c("list","sgccda", "cimdiablo")
        return(bda_object)
      } else {
        suppressMessages(mixOmics::cimDiablo(
          bda_object,
          color.blocks = paleta_niveles[1:(length(bda_object$names$blocks)-1)],
          color.Y = tr.cols[col.index.tr],
          legend.position = 0,
          size.legend = 1,
          margins = c(8, 8))
        )
        nam <- names(bda_object)
        bda_object[[length(bda_object)+1]] <- tr.cols[col.index.tr]
        bda_object[[length(bda_object)+1]] <- paleta_niveles[1:(length(bda_object$names$blocks)-1)]
        
        names(bda_object) <- c(nam, "tr.cols","lvl.cols")
        class(bda_object) <- c("list","sgccda", "cimdiablo")
        return(bda_object)
      }
      
      
    }
    
    # Network plot ----
    if (plottype == "Network") {
      if (useannot == FALSE) {
        set.seed(5881)
        
        what <- tryCatch(suppressMessages(mixOmics::network(bda_object, blocks = c(1:(length(bda_object$names$blocks)-1)), cutoff = cutoff)), error=function(e) {
          while(dev.cur()>1){dev.off()}
          texts_wizard("\nThe obtained network is too large for be displayed in RStudio plot viewer")
          texts_wizard("\nBuuuuut, it's ok, we have exported network.gml file compatible with cytoscape, and network.pdf in your working directory instead")
          #yeah! we are that good & kind
          #you are welcome
          grDevices::pdf("bda_network.pdf")
          suppressMessages(mixOmics::network(bda_object, cutoff = cutoff, blocks = c(1:(length(bda_object$names$blocks)-1))))
          dev.off() 
          grDevices::pdf("bda_network.pdf")
          a <- suppressMessages(mixOmics::network(bda_object,  cutoff = cutoff, blocks = c(1:(length(bda_object$names$blocks)-1))))
          dev.off()
          igraph::write.graph(a$gR, file = "bda_network.gml", format = "gml")
          #class(bda_object) <- c(class(bda_object), "networkplot")
          
        }
        
        )
        
        return(bda_object)
      }
      if (useannot == TRUE) {
        
        
        what <- tryCatch(suppressMessages(mixOmics::network(bda_object, blocks = c(1:length(names(bda_object))), cutoff = cutoff)), error=function(e) {
          while(dev.cur()>1){dev.off()}
          texts_wizard("\nThe obtained network is too large for be displayed in RStudio plots")
          texts_wizard("\nBuuuuut, it's ok, we have exported network.gml file compatible with cytoscape, and  network.pdf in your working directory instead")
          #yeah! we are that good & kind
          #you are welcome
          grDevices::pdf("network.pdf")
          suppressMessages(mixOmics::network(bda_list[[5]],blocks = c(1:length(names(bda_object))),cex.node.name = 0.5,cutoff = cutoff))
          dev.off() 
          grDevices::pdf("network.pdf")
          a <-suppressMessages(mixOmics::network(bda_list[[5]],blocks = c(1:length(names(bda_object))),cex.node.name = 0.5,cutoff = cutoff))
          dev.off()
          igraph::write.graph(a$gR, file = "network.gml", format = "gml")
          class(bda_list[[5]]) <-
            c(class(bda_list[[5]]), "networkplotannot")
        }
        )
        return(bda_list[[5]])
      }
    }
    # Circos plot ----
    if (plottype == "Circos") {
      
      aux.lvls <- bda_object$names$blocks[1:(length(bda_object$names$blocks)-1)]
      col.index.lvl <- unique(as.numeric(factor(aux.lvls)))
      lvl.cols <- paleta_niveles[1:length(col.index.lvl)]
      
      
      if (useannot == FALSE) {
        set.seed(5881)
        suppressMessages(mixOmics::circosPlot(bda_object, cutoff = cutoff,color.blocks = lvl.cols))
        nam <- names(bda_object)
        bda_object[[length(bda_object)+1]] <- lvl.cols
        names(bda_object) <- c(nam, "lvl.cols")
        class(bda_object) <- c( "sgccda"  , "circosplot")
        return(bda_object)
      }
      if (useannot == TRUE) {
        set.seed(5881)
        suppressMessages(mixOmics::circosPlot(bda_list[[5]], cutoff = cutoff,color.blocks = lvl.cols))
        nam <- names(bda_object)
        bda_object[[length(bda_object)+1]] <- lvl.cols
        names(bda_object) <- c(nam, "lvl.cols")
        class(bda_list[[5]]) <-
          c("sgccda" , "circosplot")
        return(bda_list[[5]])
      }
    }
  }
