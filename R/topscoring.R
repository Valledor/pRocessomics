#'@importFrom stats median
#'@importFrom utils stack

#ordenar toplot por vabs
topscoring <-
  function(object_list,
           originaldata,
           n = 20,
           mode = "abs",
           treatmentstable,
           numberofcomponents = NULL,
           treatment,
           annots) {
    options(warn = -1)
    numberofcomponents <- as.numeric(as.character(numberofcomponents))
    n <- as.numeric(as.character(n))
    if(is.list(originaldata)) {originaldata <- as.data.frame(do.call(cbind, originaldata))}
    
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
      #print("Custom palette values defined by user")
    }
    
    
    
    
  # Loadings ----
    # For PCA ----
    
      if(inherits(object_list,"pcaanalysis")){
      loads <- c()
      aux <-
        object_list$pca$rotation[, numberofcomponents] #the subset has to be like that to keep rownames
      loads <- c(loads, aux)
      a.name <- "PCA"
    }
    
    # For DA ----
      if(inherits(object_list,"daanalysis")){
      loads <- c()
      for (i in numberofcomponents) {
        aux <-
          object_list$da$loadings$X[, i] #the subset has to be like that to keep rownames
        loads <- c(loads, aux)
      }
      a.name <- "sPLS-DA"
    }
    
    # For ICA ----
      if(inherits(object_list,"icaanalysis")){
      loads <- c()
      object_list$ica$W <- t(object_list$ica$W)
      for (i in numberofcomponents) {
        aux <-
          object_list$ica$W[, i] #the subset has to be like that to keep rownames
        loads <- c(loads, aux)
      }
      a.name <- "ICA"
    }
  
    
    # For SPLS ----
      if(inherits(object_list,"splsanalysis")){
      loads <- c()
      for (i in numberofcomponents) {
        aux_x <-
          object_list$splso$loadings$X[, i] #the subset has to be like that to keep rownames
        aux_y <- object_list$splso$loadings$Y[, i]
        loads <- c(loads, aux_x, aux_y)
      }
      a.name <- "sPLS"
    }
    
    # For MCIA ----
      if(inherits(object_list,"mciaanalysis")){
      loads <- c()
      for (i in numberofcomponents) {
        aux <-
          object_list$mcia$mcoa$Tco[, i] #the subset has to be like that to keep rownames
        names(aux) <-  rownames(object_list$mcia$mcoa$Tco)
        loads <- c(loads, aux)
      }
      a.name <- "MCIA"
    }
    
    # For DIABLO ----
      if(inherits(object_list,"diabloanalysis")){
      loads <- c()
      for (j in 1:(length(object_list$diablo$loadings) - 1)) {
        loads_aux <- c()
        for (i in numberofcomponents) {
          aux <-
            (object_list$diablo$loadings[[j]])[, i] #the subset has to be like that to keep rownames
          names(aux) <- rownames(object_list$diablo$loadings[[j]])
          loads_aux <- c(loads_aux, aux)
        }
        loads <- c(loads, loads_aux)
      }
      a.name <- "DIABLO"
    }
    
  # loads sorting and so on ----
    
      loads.sort.low <- as.data.frame(sort(loads, decreasing = FALSE))
      loads.sort.high <- as.data.frame(sort(loads, decreasing = TRUE))
      loads <- as.data.frame(loads)
      loads.sort.abs <- loads[order(-abs(loads$loads)), , drop = FALSE]
     
    
    
  # subseting the top n ----
    
    if (mode == "abs") {
      topn <- as.data.frame(loads.sort.abs[1:n, ])
      rownames(topn) <- rownames(loads.sort.abs)[1:n]
      colnames(topn) <- c("Top loads abs")
      
    }
    
    if (mode == "posneg") {
      topn.high <- as.data.frame(loads.sort.high[1:n, ])
      rownames(topn.high) <- rownames(loads.sort.high)[1:n]
      colnames(topn.high) <- c("Top loads")
      topn.low <- as.data.frame(loads.sort.low[1:n, ])
      rownames(topn.low) <- rownames(loads.sort.low)[1:n]
      colnames(topn.low) <- c("Top loads")
      topn <- rbind(topn.high, topn.low)
      colnames(topn) <- c("Top loads pos&neg")
      
    }
  # get the treatment for which the median is the greatest for each topscoring variable ----
    if(treatment==4) treatment=2
    splitby <-
      factor(treatmentstable[, treatment], levels = unique(treatmentstable[, treatment]))
    spltvctr <- unique(splitby)
    
    meanvalues <-
      as.data.frame(stats::aggregate(
        originaldata[, rownames(topn)],
        by = list(splitby),
        FUN = stats::median
      ))
    meanmax <-
      do.call(rbind, apply(meanvalues, 2, function(x)
        which.max(x)))
    meanmax2 <- as.data.frame(meanvalues[meanmax, 1])
    rownames(meanmax2) <- rownames(meanmax)
    toplot <- merge(x = topn, y = meanmax2, by = "row.names")
    
  # annotations ----
    if (!is.null(annots)) {
      toplot <- merge(
        x = toplot,
        y = annots,
        by.x = "Row.names",
        by.y = "IDENTIFIER"
      )
    }
  # plot ----
    toplot <- toplot[order(-abs(toplot[, 2])), , drop = FALSE]
    if (is.null(annots)) {
      toplot[, 1] <-
        factor(toplot[, 1], levels = unique(toplot[, 1])) #to get an ordered graph
      p <-
        plotly::plot_ly(
          x = toplot[, 2],
          y = toplot[, 1],
          type = "bar",
          orientation = "h",
          color = toplot[, 3],
          colors = paleta_tratamientos[1:length(unique(spltvctr))]
        )
      p <- plotly::layout(
        p,
        xaxis = list(
          title = list(
            text = "Loadings",
            font = list(size = fontsizes[1])
          ),
          tickfont = list(size = fontsizes[2])
        ),
        
        title = list(text = (paste(a.name, "Component", numberofcomponents, "top scoring variables", sep = " ")),
                     font = list(size = fontsizes[3])),
        legend = list(font = list(size = fontsizes[4]))
      )
      
      return(p)
    }
    if (!is.null(annots)) {
      
      toplot[, 1] <-
        factor(toplot[, 1], levels = unique(toplot[, 1])) #to get an ordered graph
      p <-
        plotly::plot_ly(
          x = toplot[, 2],
          y = toplot[, 1],
          type = "bar",
          orientation = "h",
          color = toplot[, 3],
          colors = paleta_tratamientos[1:length(unique(spltvctr))],
          text = c(paste(toplot$DESCRIPTION, toplot$NAME, collapse = " ")),
          hoverinfo = 'text'
           )
      p <- plotly::layout(
        p,
        xaxis = list(
          title = list(
            text = "Loadings",
            font = list(size = fontsizes[1])
          ),
          tickfont = list(size = fontsizes[2])
        ),
        
        title = list(text = (paste(a.name, "Component", numberofcomponents, "top scoring variables", sep = " ")),
                     font = list(size = fontsizes[3])),
        legend = list(font = list(size = fontsizes[4]))
      )
      
      
      return(p)
    }
    
  }
