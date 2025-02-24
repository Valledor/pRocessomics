#' @importFrom stats kmeans

kmeans_helper<-function(matrix,i){
  set.seed(5881)
  kmeans.result <- stats::kmeans(matrix, i,iter.max=30) # Hacemos un k.means clustering (see below)
  groups <- kmeans.result$cluster #sacamos a que grupo pertenece cada variable
  withinss <- sum(kmeans.result$withinss)
  datosggplot <- as.data.frame(cbind(matrix, groups))
  datosggplot <- cbind(datosggplot, rownames(matrix))
  colnames(datosggplot)[ncol(datosggplot)] <- "ID"
  result <- list(groups,withinss,datosggplot)
  return(result)
} 
