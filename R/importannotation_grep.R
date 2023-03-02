importannotation_grep <- function(annotation){
  bincodeaux <- seq(from = 3, to = ncol(annotation), by = 2)
  nameaux <- seq(from = 4, to = ncol(annotation), by =2)
  
  # Check annotations.
  if(grepl(".", annotation[1,3], fixed = TRUE)==T){
    for(i in bincodeaux){
      annotation[, i] <- as.numeric(do.call(rbind, lapply(strsplit(as.vector(annotation[, i]), "\\."), function(x) x[[1]])))
    }
    
    for(i in nameaux){
      annotation[, i] <- as.character(do.call(rbind, lapply(strsplit(as.vector(annotation[, i]), "\\."), function(x) x[[1]])))
    }
  }
return(annotation)
  }
