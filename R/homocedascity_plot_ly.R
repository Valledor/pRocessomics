#
# Copyright 2019 Luis Valledor / Laura Lamelas
# Last revision 4.11.2019
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

#' @importFrom RColorBrewer brewer.pal
#' @importFrom plotly plot_ly layout add_trace add_annotations

homocedascity_plot_ly<-function(levenes,dataname){
  #suppressPackageStartupMessages(library("plotly"))
  #suppressPackageStartupMessages(library("RColorBrewer"))
  
  #Usar con ks.list
  lata<-lapply(levenes, function(x) length(which(x>=0.05)))
  lata<-do.call(rbind,lata)
  lata.no<-lapply(levenes, function(x) length(which(x<0.05)))
  lata.no<-do.call(rbind,lata.no)
  tata<-as.data.frame(cbind(lata,lata.no))
  colnames(tata)<-c("homocedastic","heterocedastic")
  tata<-apply(tata,1,function(x) x/sum(x)*100)
  tata<-as.data.frame(t(round(tata,2)))
  tata <- tata[nrow(tata):1,]
  #View(tata)
  colcol<-rev(RColorBrewer::brewer.pal(2,"RdYlGn"))
  plot.title.h<-paste("Homocedasticity check according to Levene test",dataname)
  r <- plotly::plot_ly()
  r <- plotly::plot_ly(tata, x = ~homocedastic, y = (factor(rownames(tata),levels = unique(rownames(tata)))), type = 'bar', orientation = 'h',
                       marker = list(color = colcol[1],
                                     line = list(color = 'rgb(248, 248, 249)', width = 1))) 
  
  r <- plotly::layout(p=r,title=plot.title.h,
                      xaxis = list(title = "",
                                   showgrid = FALSE,
                                   showline = FALSE,
                                   showticklabels = FALSE,
                                   zeroline = FALSE,
                                   domain = c(0, 100)),
                      yaxis = list(title = "",
                                   showgrid = FALSE,
                                   showline = FALSE,
                                   showticklabels = FALSE,
                                   zeroline = FALSE),
                      barmode = 'stack',
                      
                      margin = list(l = 120, r = 10, t = 140, b = 80),
                      showlegend = F)
  
  for(i in 2:ncol(tata)){
    r <- plotly::add_trace(r,x=tata[,i], y = (factor(rownames(tata),
                                                     levels = unique(rownames(tata)))),marker=list(color=colcol[[i]]))
  }
  
  r<-plotly::add_annotations(r,xref = 'paper', yref = (factor(rownames(tata),
                                                              levels = unique(rownames(tata)))), x = -0.05, y = (factor(rownames(tata),
                                                                                                                        levels = unique(rownames(tata)))),
                             xanchor = 'right',
                             text = (factor(rownames(tata),
                                            levels = unique(rownames(tata)))),
                             font = list(family = 'Arial', size = 12,
                                         color = "black"),
                             showarrow = FALSE, align = 'right')
  
  
  print(r)
}


