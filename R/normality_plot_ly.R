#Function to plot normality scores/pvalues/ahora se meter loops en plotly

#' @importFrom RColorBrewer brewer.pal
#' @importFrom plotly plot_ly layout add_trace add_annotations

normality_plot_ly<-function(ks.list,dataname){
#suppressPackageStartupMessages(library("plotly"))
#suppressPackageStartupMessages(library("RColorBrewer"))

#Usar con ks.list
count.list<-lapply(ks.list, function(x) apply(x,1,function(y) length(which(y>0.05))))
zata<-c()
for(i in 1:(ncol(ks.list[[1]])+1)){
  a<-do.call(cbind,lapply(count.list, function(x) length(which(x==(i-1)))))
  zata<-rbind(zata,a)
}
data<-as.data.frame(t(zata))  
data <- data[nrow(data):1,]
colnames(data)<-paste("x",c(1:(ncol(ks.list[[1]])+1)),sep="")

data2<-t(apply(data,1, function(x) x/sum(x)*100))
data2<-as.data.frame(data2)
colcolcol<-RColorBrewer::brewer.pal(ncol(data2),"RdYlGn")
plot.title.n<-paste("Normality check according to Kolmogorov-Smirnov test for",dataname, "dataset")
q <- plotly::plot_ly()
q <- plotly::plot_ly(data2, x = ~x1, y = (factor(rownames(data2),levels = unique(rownames(data2)))), type = 'bar', orientation = 'h',
             marker = list(color = colcolcol[1],
                           line = list(color = 'rgb(248, 248, 249)', width = 1))) 

 q<- plotly::layout(p=q,title=plot.title.n,
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

for(i in 2:ncol(data2)){
  q <-plotly::add_trace(q,x=data2[,i], y = (factor(rownames(data2),
       levels = unique(rownames(data2)))),marker=list(color=colcolcol[[i]]))
}

q<-plotly::add_annotations(q,xref = 'paper', yref =(factor(rownames(data2),
                                                   levels = unique(rownames(data2)))) , x = -0.05, y = (factor(rownames(data2),
                                                                                                               levels = unique(rownames(data2)))),
                xanchor = 'right',
                text = (factor(rownames(data2),
                               levels = unique(rownames(data2)))),
                font = list(family = 'Arial', size = 12,
                            color = "black"),
                showarrow = FALSE, align = 'right')


print(q)
}

