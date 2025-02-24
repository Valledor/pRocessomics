

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


#' @importFrom reshape2 melt 
#' @importFrom magrittr "%>%"
#' @importFrom ggplot2 ggplot aes facet_wrap geom_line stat_summary theme theme_minimal element_blank element_text
#' @importFrom plotly ggplotly layout
#' @importFrom utils txtProgressBar
#' @importFrom pheatmap pheatmap
#' @importFrom stats qt
#' @importFrom dplyr group_by summarise n mutate
#' @export
#' @name kmeans_plot
#' @title kmeans plot
#' @description A function to plot kmeans clustered groups
#' @usage kmeans_plot(kmeansobject, clusternumber = NULL, treatment = 1,
#' plottype = "kmeansclust", collapsed = FALSE,  fontsizes = c(14,10,16,12))
#' @param kmeansobject List containing kmeans analysis result generated with pRocessomics::kmeans_analysis()
#' @param clusternumber Numeric, Number of clusters to be represented. Note that the number of clusters has to be concordant with the respective kmeansobject
#' If set to NULL or not declared and selected plottype is kmeansclust all calculated plot will be displayed and in the case of heatmap plottype all variables will be displayed as heatmap rows.
#' @param treatment Numeric value to determine the treatments to be displayed in the x axes. Available options are:
  #' \itemize{
  #' \item 1 according to treatment 1 (treatment1col)
  #' \item 2 according to treatment 2 (treatment2col)
  #' \item 3 according to the combination of both treatments
  #' \item 4 multiple-kmeans: X-axes according to treatment 2 and different lines for treatment 1
  #' }
#' @param plottype Character, class of the plot to display available options are kmeansclust and heatmap
#' @param collapsed Boolean, TRUE indicates multiple overlap depiction
#' @param fontsizes Vector (length 4) containing the font sizes to be used in the plot
#' @return kmeans plot
#' @author Luis Valledor and Laura Lamelas
#'
kmeans_plot <-
  function(kmeansobject,
           clusternumber = NULL,
           treatment = 1,
           plottype = "kmeansclust",
           collapsed = FALSE,
           fontsizes = c(14, 10, 16, 12)) {
    
# Flags, checks, aesthetics ----
    variable <- NULL
    value <- NULL
    ID <- NULL
    groups <- NULL
    CI_lower <- CI_upper <- sem <- NULL
    
    
    if (is.null(treatment)) treatment <- 1
    options(warn = -1)
   if(treatment == 4){
    treatmentstable <- unique(kmeansobject$treatments)
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
    AAcolor = factor(treatmentstable[, 1], levels = unique(treatmentstable[, 1]))
    AAcolors = paleta_tratamientos[1:length(unique(treatmentstable[, 1]))]
    
   }}
    
   
    if(!inherits(kmeansobject, "kmeansanalysis"))  
      stop("A kmeansanalysis object is required. Please run kmeans_analysis first")
# kmeans facets plot ----
    if(plottype == "kmeansclust"){
      # treatment %in% 1:3 singlekmeans ----
    if(treatment %in% 1:3){
        if (!is.null(clusternumber)) {
        clusteringsdone <-
          as.numeric(gsub(" Clusters", "", names(kmeansobject$kmeans_list), fixed = TRUE)) #sacamos los numeros de los clusters hechos
        if (clusternumber %in% clusteringsdone == FALSE)
          stop(
            "The used kmeans analysis was performed with different numbers of clusters, please re-run kmeans_analysis() function"
          ) #vemos que ese numero de groups se haya hecho de forma efectiva
        
        cluster <-
          which(as.numeric(gsub(
            " Clusters", "", names(kmeansobject$kmeans_list), fixed = TRUE
          )) == clusternumber) #sacamos que posicion en la lista de kmeans list ocupa el numero de clusters deseado
        datosggplot <-
          kmeansobject$kmeans_list[[cluster]] #sacamos los datos requeridos
        datosggplot <-
          reshape2::melt(datosggplot, (ncol(datosggplot) - 1):(ncol(datosggplot))) #melt para ponerlo en formato para ggplot2
        
        datosggplot <- datosggplot[order(datosggplot$groups), ]
        datosggplot2 <- datosggplot
        datosggplot2[, 1] <- paste("Cluster", datosggplot2[, 1])
        a <-
          factor(datosggplot2[, 1], levels = c(paste(
            "Cluster", unique(datosggplot[, 1]), sep = " "
          )))
        myplot <-
          ggplot2::ggplot(
            datosggplot2,
            ggplot2::aes(
              x = variable,
              y = value,
              group = ID,
              colour = as.factor(groups),
              text = datosggplot2$ID
            )
          ) +
          ggplot2::facet_wrap(facets = a , scales = "free") +
          ggplot2::geom_line(linetype = "dotted",
                             alpha = 1 / 5,
                             lwd = 1) +
          ggplot2::stat_summary(aes(group = groups), fun.y = mean, geom = "line") +
          ggplot2::theme_minimal() + ggplot2::theme(
            legend.position = "none",
            axis.title.x = ggplot2::element_blank(),
            axis.title.y = element_blank(),
            axis.text.x = ggplot2::element_text(
              angle = 45,
              vjust = 1,
              hjust = 1
            )
          )
        
        gg <- plotly::ggplotly(myplot, tooltip = c("text"))
        
        gg <-
          plotly::layout(
            gg,
            title = "K-means grouping",
            titlefont = list(size = fontsizes[3]),
            legend = list(font = list(size = fontsizes[4]))
          )
        #class(gg)<-c("plotly","htmlwidget","singlekmeansplot")
        print(gg)
        singlekmeansplot <- list(gg, myplot)
        class(singlekmeansplot) <- "singlekmeansplot"
        return(singlekmeansplot)
      }
      ## generate plots for all the clusters 
      if (is.null(clusternumber)){
        ggplotlist <- c()
        myplotlist <- c()
        cat("Generating k-means plots. It may take a while\n")
        pb <-
          utils::txtProgressBar(
            min = 0,
            max = length(kmeansobject$kmeans_list),
            style = 3
          )
        for (i in 1:length(kmeansobject$kmeans_list)) {
          datosggplot <-
            kmeansobject$kmeans_list[[i]] #sacamos los datos requeridos
          datosggplot <-
            reshape2::melt(datosggplot, (ncol(datosggplot) - 1):(ncol(datosggplot))) #melt para ponerlo en formato para ggplot2
         
          datosggplot <- datosggplot[order(datosggplot$groups), ]
          datosggplot2 <- datosggplot
          datosggplot2[, 1] <- paste("Cluster", datosggplot2[, 1])
          a <-
            factor(datosggplot2[, 1], levels = c(paste(
              "Cluster", unique(datosggplot[, 1]), sep = " "
            )))
          myplot <-
            ggplot2::ggplot(
              datosggplot2,
              aes(
                x = variable,
                y = value,
                group = ID,
                colour = as.factor(groups),
                text = datosggplot2$ID
              )
            ) +
            facet_wrap(facets = a , scales = "free") +
            geom_line(linetype = "dotted",
                      alpha = 1 / 5,
                      lwd = 1) +
            stat_summary(aes(group = groups), fun.y = mean, geom = "line") +
            theme_minimal() + theme(
              legend.position = "none",
              axis.title.x = element_blank(),
              axis.title.y = element_blank(),
              axis.text.x = element_text(
                angle = 45,
                vjust = 1,
                hjust = 1
              )
            )
          
          gg <- plotly::ggplotly(myplot, tooltip = c("text"))
          
          gg <-
            plotly::layout(
              gg,
              title = "K-mean grouping",
              titlefont = list(size = fontsizes[3]),
              legend = list(font = list(size = fontsizes[4]))
            )
          print(gg)
          myplotlist[[i]] <-
            myplot #esto es una GRAN chapuza, pero no se exportar de otra forma
          ggplotlist[[i]] <- gg
          utils::setTxtProgressBar(pb, i)
          
        }
        names(ggplotlist) <- names(kmeansobject$kmeans_list)
        names(myplotlist) <- names(kmeansobject$kmeans_list)
        #class(ggplotlist)<-c("plotly","htmlwidget","multiplekmeansplot")
        kmeansmultiple <- list(ggplotlist, myplotlist)
        class(kmeansmultiple) <- "kmeansmultipleplot"
        return(kmeansmultiple)
      }
      }
      # treatment==4 Overlapped kmeans ----
    if(treatment == 4){ 
       if (!is.null(clusternumber)) {
        clusteringsdone <-
          as.numeric(gsub(" Clusters", "", names(kmeansobject$kmeans_list), fixed = TRUE)) #sacamos los numeros de los clusters hechos
        if (clusternumber %in% clusteringsdone == FALSE)
          stop(
            "The used kmeans analysis was performed with different numbers of clusters, please re-run kmeans_analysis() function"
          ) #vemos que ese numero de groups se haya hecho de forma efectiva
        
        cluster <-
          which(as.numeric(gsub(
            " Clusters", "", names(kmeansobject$kmeans_list), fixed = TRUE
          )) == clusternumber) #sacamos que posicion en la lista de kmeans list ocupa el numero de clusters deseado
        datosggplot <-
          kmeansobject$kmeans_list[[cluster]] #sacamos los datos requeridos
        datosggplot <-
          reshape2::melt(datosggplot, (ncol(datosggplot) - 1):(ncol(datosggplot))) #melt para ponerlo en formato para ggplot2
      
         datosggplot <- datosggplot[order(datosggplot$groups), ]
        datosggplot2 <- datosggplot
        datosggplot2[, 1] <- paste("Cluster", datosggplot2[, 1])
        a <-
          factor(datosggplot2[, 1], levels = c(paste(
            "Cluster", unique(datosggplot[, 1]), sep = " "
          )))
        ######caution!-----
        key <- unique(kmeansobject$treatments) #required treatment combinations
        
        ## replace column variable with the corresponding columns of key
        row.ind <- match(datosggplot2$variable, key$V3)
        t1 <- key$V1[row.ind]
        t2 <- key$V2[row.ind]
        datosggplot3 <- cbind(datosggplot2[,1:2], "t1"= t1, "t2"=t2,"value"=datosggplot2[,4] )
        
        ########  ribbon-----
        if (collapsed == T){
          #sosorted <- unique(datosggplot3$groups) 
          df_tidy_mean <- datosggplot3 %>%
            dplyr::group_by(groups, t1, t2) %>%
            dplyr::summarise(n = dplyr::n(),
                      mean = mean(value),
                   
                      sd = sd(value),
                   ) %>%
            dplyr::mutate(sem = sd / sqrt(n),
                   CI_lower = mean + stats::qt((1-0.95)/2, n) * sem,
                   CI_upper = mean - stats::qt((1-0.95)/2, n) * sem)
          
          dff2 <- as.data.frame(df_tidy_mean)
          dff <- dff2[order(unlist(sapply(dff2$groups, function(x) which(unique(datosggplot3$groups) ==x)))),]
          
          dff$groups <- paste(dff$groups," (",dff$n,")", sep="")
          dff$t2 <- factor(dff$t2, levels = unique(t2))
          dff_sorted <- dff %>%
            dplyr::arrange(groups, t1, t2)
          
            myNEWplot <-
            ggplot2::ggplot(
              dff,
              ggplot2::aes(
                x = t2,
                y = mean,
                group = t1,
                color = t1
                ),
              text = t1
            ) +
            
            ggplot2::scale_fill_manual(values = AAcolors) + #opacity
            ggplot2::scale_color_manual(values = AAcolors) + #standard
            ggplot2::facet_wrap(facets = factor(dff$groups, levels= unique(dff$groups)) , scales = "free") +
            ggplot2::geom_ribbon(aes(ymin=CI_lower,ymax=CI_upper,fill=t1)#, colour=t1, 
                                 ,alpha = 0.3
                                 ) +

            ggplot2::geom_line(aes(x=t2, y= dff$mean, color=t1
                                   #,
                                   #colours = AAcolors
                                   )) +
            
            ggplot2::theme_minimal() + ggplot2::theme(
              legend.position = "top",
              axis.title.x = ggplot2::element_blank(),
              axis.title.y = element_blank(),
              axis.text.x = ggplot2::element_text(
                angle = 0,
                vjust = 1,
                hjust = 1
              )
            )
          
          #myNEWplot
          
          gg <- plotly::ggplotly(myNEWplot, tooltip = c("text"))
          gg <-
            plotly::layout(
              gg,
              title = list(text="K-means clustering", y=1.5, x=0.05),
              titlefont = list(size = fontsizes[3]),
              showlegend = F,
              legend = list(orientation = "v", traceorder = "grouped")
              # yanchor ="top")
              # legend = list(font = list(size = fontsizes[4]))
            )
          
          print(gg)
          singlekmeansplot <- list(gg, myNEWplot)
          class(singlekmeansplot) <- "singlekmeansplot"
          return(singlekmeansplot)  
        }
        
        ########   sparse lines----  
        if (collapsed == F){
          myplot <-
            ggplot2::ggplot(
              datosggplot3,
              ggplot2::aes(
                x = t2,
                y = value,
                group = ID,
                colour = as.factor(datosggplot3$t1),
                text = datosggplot3$t1
              )
            ) +
            
            #ggplot2::scale_color_manual(values=as.factor(datosggplot3$t1)) +
            ggplot2::facet_wrap(facets = a , scales = "free") +
            ggplot2::geom_line(linetype = "solid",
                               alpha = 0.2,
                               lwd = 0.4) +
            ggplot2::stat_summary(aes(group = t1), fun.y = mean, geom = "line") +
            ggplot2::theme_minimal() + ggplot2::theme(
              legend.position = "none",
              axis.title.x = ggplot2::element_blank(),
              axis.title.y = element_blank(),
              axis.text.x = ggplot2::element_text(
                angle = 0,
                vjust = 1,
                hjust = 1
              )
            )
          
          gg <- plotly::ggplotly(myplot, tooltip = c("text"))
          
          gg <-
            plotly::layout(
              gg,
              legend = list(font = list(size = fontsizes[4])),
              title = list(text="Multiple K-means clustering", y=1.5, x=0.05),
              titlefont = list(size = fontsizes[3]),
              showlegend = T,
              legend = list(orientation = "h", y=1.1, x = 0.6)
            )
          
          
          print(gg)
          singlekmeansplot <- list(gg, myplot)
          class(singlekmeansplot) <- "singlekmeansplot"
          return(singlekmeansplot)  
        }  
          
        
      }
      ## generate plots for all the clusterings ----
      if (is.null(clusternumber)){
        ggplotlist <- c()
        myplotlist <- c()
        cat("Generating k-means plots. It may take a while\n")
        pb <-
          utils::txtProgressBar(
            min = 0,
            max = length(kmeansobject$kmeans_list),
            style = 3
          )
        for (i in 1:length(kmeansobject$kmeans_list)) {
         
          
           datosggplot <-
            kmeansobject$kmeans_list[[i]] #sacamos los datos requeridos
          datosggplot <-
            reshape2::melt(datosggplot, (ncol(datosggplot) - 1):(ncol(datosggplot))) #melt para ponerlo en formato para ggplot2
  
          datosggplot <- datosggplot[order(datosggplot$groups), ]
          datosggplot2 <- datosggplot
          datosggplot2[, 1] <- paste("Cluster", datosggplot2[, 1])
          a <-
            factor(datosggplot2[, 1], levels = c(paste(
              "Cluster", unique(datosggplot[, 1]), sep = " "
            )))
          #### treatment 1 to 3----        
          if(treatment %in% 1:3){
            
            myplot <-
              ggplot2::ggplot(
                datosggplot2,
                aes(
                  x = variable,
                  y = value,
                  group = ID,
                  colour = as.factor(groups),
                  text = datosggplot2$ID
                )
              ) +
              facet_wrap(facets = a , scales = "free") +
              geom_line(linetype = "dotted",
                        alpha = 1 / 5,
                        lwd = 1) +
              stat_summary(aes(group = groups), fun.y = mean, geom = "line") +
              theme_minimal() + theme(
                legend.position = "none",
                axis.title.x = element_blank(),
                axis.title.y = element_blank(),
                axis.text.x = element_text(
                  angle = 45,
                  vjust = 1,
                  hjust = 1
                )
              )
            
            gg <- plotly::ggplotly(myplot, tooltip = c("text"))
            
            gg <-
              plotly::layout(
                gg,
                title = "Multiple k-means grouping",
                titlefont = list(size = fontsizes[3]),
                legend = list(font = list(size = fontsizes[4]))
              )
            print(gg)
          }
          
          #### treatment == 4----
          if(treatment == 4){
            
            ######caution!-----
            key <- unique(kmeansobject$treatments) #required treatment combinations
            
            ## replace column variable with the corresponding columns of key
            row.ind <- match(datosggplot2$variable, key$V3)
            t1 <- key$V1[row.ind]
            t2 <- key$V2[row.ind]
            datosggplot3 <- cbind(datosggplot2[,1:2], "t1"= t1, "t2"=t2,"value"=datosggplot2[,4] )
            
            ########  ribbon-----
            if (collapsed == T){
              #sosorted <- unique(datosggplot3$groups) 
              df_tidy_mean <- datosggplot3 %>%
                dplyr::group_by(groups, t1, t2) %>%
                dplyr::summarise(n = dplyr::n(),
                          mean = mean(value),
                          
                          sd = sd(value),
                ) %>%
                dplyr::mutate(sem = sd / sqrt(n),
                       CI_lower = mean + stats::qt((1-0.95)/2, n) * sem,
                       CI_upper = mean - stats::qt((1-0.95)/2, n) * sem)
              
              dff2 <- as.data.frame(df_tidy_mean)
              dff <- dff2[order(unlist(sapply(dff2$groups, function(x) which(unique(datosggplot3$groups) ==x)))),]
              dff$groups <- paste(dff$groups," (",dff$n,")", sep="")
           
              myplot <-
                ggplot2::ggplot(
                  dff,
                  ggplot2::aes(
                    x = t2,
                    y = mean,
                    group = t1,
                    color = t1),
                  text = t1
                ) +
                
                #ggplot2::scale_fill_manual(as.factor(datosggplot3$t1)) + #opacity
                #ggplot2::scale_color_manual(as.factor(datosggplot3$t1)) + #standard
                ggplot2::facet_wrap(facets = factor(dff$groups, levels= unique(dff$groups)) , scales = "free") +
                ggplot2::geom_ribbon(aes(ymin=CI_lower,ymax=CI_upper,fill=t1)#, colour=t1, 
                                     ,alpha = 0.3
                ) +
                
                ggplot2::geom_line(aes(x=t2, y= dff$mean, color=t1)) +
                
                ggplot2::theme_minimal() + ggplot2::theme(
                  legend.position = "top",
                  axis.title.x = ggplot2::element_blank(),
                  axis.title.y = element_blank(),
                  axis.text.x = ggplot2::element_text(
                    angle = 0,
                    vjust = 1,
                    hjust = 1
                  )
                )
              
              #myNEWplot
              
              gg <- plotly::ggplotly(myplot, tooltip = c("text"))
              gg <-
                plotly::layout(
                  gg,
                  title = list(text="K-means clustering", y=1.5, x=0.05),
                  titlefont = list(size = fontsizes[3]),
                  showlegend = T,
                  legend = list(orientation = "h", y=1.1, x = 0.6)
                )
              
             
            }
            
            ########   sparse lines----  
            if (collapsed == F){
              myplot <-
                ggplot2::ggplot(
                  datosggplot3,
                  ggplot2::aes(
                    x = t2,
                    y = value,
                    group = ID,
                    colour = as.factor(datosggplot3$t1),
                    text = datosggplot3$t1
                  )
                ) +
                ggplot2::facet_wrap(facets = a , scales = "free") +
                ggplot2::geom_line(linetype = "solid",
                                   alpha = 0.2,
                                   lwd = 0.4) +
                ggplot2::stat_summary(aes(group = t1), fun.y = mean, geom = "line") +
                ggplot2::theme_minimal() + ggplot2::theme(
                  legend.position = "none",
                  axis.title.x = ggplot2::element_blank(),
                  axis.title.y = element_blank(),
                  axis.text.x = ggplot2::element_text(
                    angle = 0,
                    vjust = 1,
                    hjust = 1
                  )
                )
              
              gg <- plotly::ggplotly(myplot, tooltip = c("text"))
              
              gg <-
                plotly::layout(
                  gg,
                  legend = list(font = list(size = fontsizes[4])),
                  title = list(text="Multiple K-means clustering", y=1.5, x=0.05),
                  titlefont = list(size = fontsizes[3]),
                  showlegend = T,
                  legend = list(orientation = "h", y=1.1, x = 0.6)
                )
              
              
             
            }  
            
          }  
          myplotlist[[i]] <-
            myplot #esto es una GRAN chapuza, pero no se exportar de otra forma
          ggplotlist[[i]] <- gg
          utils::setTxtProgressBar(pb, i)
          
        }
        
        
        names(ggplotlist) <- names(kmeansobject$kmeans_list)
        names(myplotlist) <- names(kmeansobject$kmeans_list)
        #class(ggplotlist)<-c("plotly","htmlwidget","multiplekmeansplot")
        kmeansmultiple <- list(ggplotlist, myplotlist)
        class(kmeansmultiple) <- "kmeansmultipleplot"
        return(kmeansmultiple)
      }
      }
    }
    
   
    # heatmap ----
    if (plottype == "heatmap") {
      if(is.null(clusternumber)) clusternumber <- nrow(unique(t(kmeansobject$kmeans_matrix)))
      #print(clusternumber)
      values <- kmeansobject$kmeans_matrix
      set.seed(151186) #required to reproducibility
      p <-
        pheatmap::pheatmap(
          t(values),
          color = paleta_continuo,
          kmeans_k = clusternumber,
          annotation_names_col = F,
          annotation_names_row = F,
          annotation_row = NULL,
          annotation_col = NULL,
          show_rownames = T,
          #clustering_distance_rows = hmparameters[1],
          #clustering_distance_cols = hmparameters[2],
          #clustering_method = hmparameters[3],
          display_numbers = FALSE#,
          #scale = hmparameters[5]
        )
      class(p) <- "Plot_kmeans_heatmap"
      #return(p) #ojo, para plotearlo a antojo posteriormente usar grid::grid.draw(p$gtable)
    }
  }
