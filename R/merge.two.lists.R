merge.two.lists<-function(list1,list2){
  aux.names <- unique(c(names(list1),names(list2)))
  to.combine <- names(list2[which(names(list2) %in% names(list1))])
  to.append.list1 <- list1[(names(list1)[!names(list1) %in% names(list2)])]
  to.append.list2 <- list2[(names(list2)[!names(list2) %in% names(list1)])]
  common.list1 <- list1[(to.combine)]
  common.list2 <- list2[(to.combine)]
  common.list <-list()
  for (i in 1:length(to.combine)){
    common.list[[i]] <- (as.vector(c(common.list1[[i]],common.list2[[i]])))
  }
  names(common.list) <- to.combine
  final.list <- c(common.list,to.append.list1,to.append.list2)
}
