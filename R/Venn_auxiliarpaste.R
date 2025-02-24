venn_auxiliarpaste<-function(v){  
  v<-v[!is.na(v)]
  v<-paste(v,collapse="_||_")
  return(v)
}
