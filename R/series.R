series <- function(x){
  x<-ncol(x)
  if(x>300) if(x>250) y<- c(seq(10,50,5),seq(75,150,25),200,250,300)
  if(x<300) if(x>250) y<- c(seq(10,50,5),seq(75,150,25),200,250,x)
  if(x>=100&x<=250) y<-c(seq(10,50,5),round(seq(75,x,10)))
  if(x>=25&x<100) y<-unique(c(round(seq(10,x,13))))      
  if(x<25) y <- unique(c(round(seq(2,x,5))))
  return(y)
}
