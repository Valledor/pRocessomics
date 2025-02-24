#Function to plot normality scores/pvalues/ahora se meter loops en plotly----------

#' @importFrom grDevices colorRampPalette
#' @importFrom RColorBrewer brewer.pal
#' @importFrom plotly plot_ly layout add_trace add_annotations subplot

#normality
parametricity_plot_ly <- function(ks.list, levenes, dataname) {
  
  #Normalitity check
  count.list <-
    lapply(ks.list, function(x)
      apply(x, 1, function(y)
        length(which(y > 0.05))))
  zata <- c()
  for (i in 1:(ncol(ks.list[[1]]) + 1)) {
    a <-
      do.call(cbind, lapply(count.list, function(x)
        length(which(x == (
          i - 1
        )))))
    zata <- rbind(zata, a)
  }
  data <- as.data.frame(t(zata))
  data <- data[nrow(data):1, ]
  colnames(data) <- paste("x", c(1:(ncol(ks.list[[1]]) + 1)), sep = "")
  
  data2 <- t(apply(data, 1, function(x)
    x / sum(x) * 100))
  data2 <- as.data.frame(data2)
  colcolcol <-
    grDevices::colorRampPalette(RColorBrewer::brewer.pal(ncol(data2), "RdYlGn"))(ncol(data2))
  plot.title.n <-
    paste("Normality check according to Kolmogorov-Smirnov test for",
          dataname,
          "dataset")
  
  q <-
    plotly::plot_ly(
      data2,
      x = ~ x1,
      y = (factor(rownames(data2), levels = unique(rownames(
        data2
      )))),
      type = 'bar',
      orientation = 'h',
      name = "none",
      marker = list(
        color = colcolcol[1],
        lengendgroup = "normality",
        line = list(
          color = 'rgb(248, 248, 249)',
          width = 1,
          name = "none"
        )
      )
    )
  
  for (i in 2:ncol(data2)) {
    trace.names = paste(i - 1, "treats")
    trace.names[i == ncol(data2)] <- "all"
    q <-
      plotly::add_trace(
        q,
        x = data2[, i],
        name = trace.names,
        lengendgroup = "normality",
        y = (factor(rownames(data2),
                    levels = unique(
                      rownames(data2)
                    ))),
        marker = list(color = colcolcol[[i]])
      )
    q <-  plotly::layout(
      q,
      title = plot.title.n,
      xaxis = list(
        title = "",
        showgrid = FALSE,
        showline = FALSE,
        showticklabels = FALSE,
        zeroline = FALSE,
        domain = c(0, 100)
      ),
      yaxis = list(
        title = "",
        showgrid = FALSE,
        showline = FALSE,
        showticklabels = FALSE,
        zeroline = FALSE
      ),
      barmode = 'stack',
      margin = list(
        l = 120,
        r = 10,
        t = 140,
        b = 80
      ),
      legend = list(
        orientation = 'h',
        x = 0,
        y = 4,
        traceorder = "normal",
        xanchor = 'left',
        tracegroupgap = 12
      )
    )
    
    
  }
  
  q <-
    plotly::add_annotations(
      q,
      xref = 'paper',
      yref = (factor(rownames(data2), levels = unique(rownames(
        data2
      )))) ,
      x = -0.05,
      y = (factor(rownames(data2), levels = unique(rownames(
        data2
      )))),
      xanchor = 'right',
      text = (factor(rownames(data2), levels = unique(rownames(
        data2
      )))),
      font = list(
        family = 'Arial',
        size = 12,
        color = "black"
      ),
      showarrow = FALSE,
      align = 'right'
    )
  
  
  #Homocedasticy check
  lata <- lapply(levenes, function(x)
    length(which(x >= 0.05)))
  lata <- do.call(rbind, lata)
  lata.no <- lapply(levenes, function(x)
    length(which(x < 0.05)))
  lata.no <- do.call(rbind, lata.no)
  tata <- as.data.frame(cbind(lata, lata.no))
  colnames(tata) <- c("homocedastic", "heterocedastic")
  tata <- apply(tata, 1, function(x)
    x / sum(x) * 100)
  tata <- as.data.frame(t(round(tata, 2)))
  tata <- tata[nrow(tata):1, ]
  colcol <- rev(RColorBrewer::brewer.pal(2, "RdYlGn"))
  plot.title.h <-
    paste("Homocedasticity check according to Levene test",
          dataname,
          "dataset")
  r <-
    plotly::plot_ly(
      tata,
      x = ~ homocedastic,
      lengendgroup = "Variance",
      y = (factor(rownames(tata), levels = unique(rownames(
        tata
      )))),
      type = 'bar',
      orientation = 'h',
      name = "Homocedastic",
      marker = list(
        color = colcol[1],
        line = list(color = 'rgb(248, 248, 249)', width = 1)
      )
    )
  
  r <-
    plotly::add_trace(
      r,
      x = tata[, 2],
      name = "Heterocedastic",
      lengendgroup = "Variance",
      y = (factor(rownames(tata),
                  levels = unique(rownames(
                    tata
                  )))),
      marker = list(color = colcol[2])
    )
  
  r <- plotly::layout(
    r,
    title = plot.title.h,
    xaxis = list(
      title = "",
      showgrid = FALSE,
      showline = FALSE,
      showticklabels = FALSE,
      zeroline = FALSE,
      domain = c(0, 100)
    ),
    yaxis = list(
      title = "",
      showgrid = FALSE,
      showline = FALSE,
      showticklabels = FALSE,
      zeroline = FALSE
    ),
    barmode = 'stack',
    margin = list(
      l = 120,
      r = 10,
      t = 140,
      b = 80
    ),
    legend = list(
      orientation = 'h',
      x = 0,
      y = 4,
      traceorder = "normal",
      xanchor = 'left'
    )
  )
  
  
  #print(q)
  #print(r)
  #Subplot
  margen <- 0.08
  a <- nrow(tata) * 0.1 + 0.1
  b <- (a + margen) * 2
  subtitulo <- paste(dataname, "datalist")
  locations <- (seq(1:nrow(tata)) - 1) * (0.1) + 0.05
  p <-
    plotly::subplot(
      q,
      r,
      nrows = 2,
      shareX = T,
      which_layout = 'merge',
      margin = margen,
      titleY = T
    )
  p <-
    plotly::layout(
      p = p,
      title = paste(
        "Parametricity tests according to Kolmogorov-Smirnov and Levene\n",
        subtitulo
      ),
      annotations = list(
        list(
          x = 0,
          y = b,
          text = "Normal distribution check",
          showarrow = F,
          xref = 'paper',
          yref = 'paper'
        ),
        list(
          x = 0,
          y = a,
          text = "Homocedascity check",
          showarrow = F,
          xref = 'paper',
          yref = 'paper'
        )
      ),
      legend = list(
        x = 0,
        y = -0.1,
        legendgroup = c("normality", "variance"),
        tracegroupgap = 1
      )
    )
  
  
  
  p <-
    plotly::add_annotations(
      p,
      xref = 'paper',
      yref = 'paper',
      x = -0.05,
      y = locations,
      xanchor = 'right',
      yanchor = 'middle',
      text = (factor(rownames(tata), levels = unique(rownames(
        tata
      )))),
      font = list(
        family = 'Arial',
        size = 12,
        color = "black"
      ),
      showarrow = FALSE,
      align = 'right'
    )
  print(p)
}
