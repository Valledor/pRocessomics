sumpergroup <-function(lista,grupos){
  p1<-lapply(grupos,function(x) as.data.frame(lista[,colnames(lista)%in%x]))
  p2 <-lapply(p1, function(x) rowSums(x,na.rm=T))
  p3 <- do.call(cbind,p2)
  return(p3)
}