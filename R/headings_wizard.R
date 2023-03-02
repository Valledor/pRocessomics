headings_wizards<-function(title.text){
  screenW<-getOption("width")
  all<-(c("\n",rep("#",screenW),sep = ""))
  title.length<-nchar(title.text)
  if(screenW<title.length+2){
    heading.text<- strsplit(title.text," ")
    h.text<-as.vector(do.call(rbind, heading.text))
    cat(all,sep="")
    cat(paste0(c("\n","#",h.text,"#",sep="")),fill=T)
    cat(all,sep="")
  }
  else{
  spaces<-(screenW-title.length-2)/2
  cat(all,sep="")
  middle<-cat(c("\n","#",rep(" ",spaces),title.text,rep(" ",spaces),"#"),sep="")
  cat(all,sep="")
  }
}

texts_wizard<-function(body.text){ #DON'T USE \n
  text1 <- strsplit(body.text," ")
  text1<-as.vector(do.call(rbind, text1))
  cat(paste0(text1,sep=""),fill=T)
}
