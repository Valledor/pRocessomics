#Circos matrix for circos plot...

#### Main: Circos matrix---------------
Circos_matrix<-function(datalist, initialrow, initialcolumn, threshold){
  if(class(datalist)!="POL"){stop ("Error, please re-run preprocess steps")} #checking data is preprocessed
  n<-length(datalist)
  #let's supose threshold can be a single value or a vector of length n^2
  if(length(threshold)==1){
    t.vec<-rep(threshold,n^2)
    thres<-matrix(t.vec,nrow=n,ncol=n, byrow = TRUE)
  }
  else if(length(threshold)==n^2){
    thres<-matrix(thres,nrow=n,ncol=n, byrow = TRUE)
  }
  else {stop ("threshold should be a single number or a vector of length n^2, being n the number of omic levels to integrate") }
  ###other option (the one I dislike the most...)
  #threshold is a vector of two components, the first one for smbn and the second for spls
  #if(length(threshold)==2){
  #  thres_a<-threshold[1]
  #  thres_b<-threshold[2]
  #}
  
  links.No<-data.frame()
  
  for (i in 1:n){
    for(j in 1:n){
      if (i != j){
        
        comps<-length(initialrow:nrow(datalist[[i]]))-(ceiling(length(initialrow:nrow(datalist[[i]]))/10)+1) #maximum allowed ncomp
        
        Xn<-as.matrix(datalist[[i]][initialrow:nrow(datalist[[i]]),initialcolumn:ncol(datalist[[i]])])
        Yn<-as.matrix(datalist[[j]][initialrow:nrow(datalist[[j]]),initialcolumn:ncol(datalist[[j]])])
        spls.ij<-spls(Xn,Yn,ncomp=comps, mode="regression")
        
        set.seed(1511)
        nfolds<-length(initialrow:nrow(datalist[[i]]))
        
        A<-perf(spls.ij, validation = "Mfold", folds = nfolds, progressBar = FALSE, nrepeat = 100)
        if(all(A$Q2.total[,1]>=0.0975)){
          aux_comp<-1:comps
        }
        else if(all(A$Q2.total[,1]<0.0975)){
          next
        }
        else{
          # k=1
          #  while(A$Q2.total[K,1]>=0.0975){ #rule of thumbs of mixOmics...
          #    k=k+1   #esto seria correcto si los componentes del spls tuvieran un Q2.total creciente pero no es asi (siempre), solo en el puto ejemplo
          aux_comp<-1:min(which(A$Q2.total>=0.0975)) #max or min??????
          
        }
        
        net.result.ij<- new_network2(mat=spls.ij, comp = aux_comp, cutoff = thres[i,j])  
        mat1<-as.matrix(net.result.ij)
        mat2<-melt(mat1,id.vars=rownames(mat1),measure.vars=colnames(mat1))
        colnames(mat2)<-c("Source","Target","Weight")
        mat3<-mat2[abs(mat2$weight)>=thres[i,j],]
        links_cutoff<-nrow(mat3)
        links.No[i,j]<-links_cutoff
      }
      else{
        single.num.matrix<-datalist[[i]][initialrow:nrow(datalist[[i]]),initialcolumn:ncol(datalist[[i]])]
        smbn_i<-smbn(single.num.matrix,thres[i,j])
        autolink<-nrow(smbn_i)
        links.No[i,j]<-autolink
      }
    }
    
  }
  return (links.No)
}



#####Packages-------------
#Package requirements
Package_aux<-function(required.packages){
  #required.packages is a vector made of the quoted names of the required packages
  #auxiliary function for checking, installing and attaching packages
  i=1
  while(i <= (length(required.packages))){
    if(required.packages[i] %in% rownames(installed.packages())==FALSE){
      
      install.packages(required.packages[i])
      b<-" installed and "
    }
    else{b<-" "}
    suppressPackageStartupMessages(library(required.packages[i],character.only = TRUE))
    i=i+1
  }
  print(paste(required.packages," has been",b,"loaded.",sep=""))
}

needed.packages<-c("mixOmics","reshape2","plyr","dplyr","plotly")
Package_aux(needed.packages)

#input:
# preprocessed datalist, thresholds, initialrow initialcolumn don't need treatments here 

#### Aux_functions--------
new_network2<-function (mat, comp = NULL, cutoff = NULL,col.names=TRUE) {
  arg.call = match.call()
  user.arg = names(arg.call)[-1]
  err = tryCatch(mget(names(formals()), sys.frame(sys.nframe())), 
                 error = function(e) e)
  if ("simpleError" %in% class(err)) 
    stop(err[[1]], ".", call. = FALSE)
  function.arg = names(mget(names(formals()), sys.frame(sys.nframe())))
  not.arg = !(user.arg %in% function.arg)
  if (any(not.arg)) {
    unused.arg = user.arg[not.arg]
    not.arg = which(not.arg) + 1
    output = rep("", length(not.arg))
    for (i in 1:length(not.arg)) {
      output[i] = paste0(unused.arg[i], " = ", arg.call[[not.arg[i]]])
    }
    output = paste0("(", paste(output, collapse = ", "), 
                    ").")
    msg = "unused argument "
    if (length(not.arg) > 1) 
      msg = "unused arguments "
    stop(msg, output, call. = FALSE)
  }
  class.object = class(mat)
  object.pls = c("pls", "spls", "mlspls")
  object.rcc = "rcc"
  object.blocks = c("sgcca", "rgcca")
  if (!any(class.object %in% c(object.pls, object.rcc, object.blocks, 
                               "matrix"))) 
    stop(" 'network' is only implemented for the following objects: matrix, pls, plsda, spls, splsda, rcc, sgcca, rgcca, sgccda", 
         call. = FALSE)
  if (any(class.object %in% c(object.rcc, object.pls))) {
    p = ncol(mat$X)
    q = ncol(mat$Y)
    n = nrow(mat$X)
    ncomp = mat$ncomp
    if (is.null(comp)) 
      comp = 1:mat$ncomp
    if (length(comp) == 1) {
      if (comp > ncomp) {
        stop("the elements of 'comp' must be smaller than or equal to ", 
             ncomp, ".", call. = FALSE)
      }
      else if (!is.numeric(comp) || comp <= 0) {
        stop("invalid value for 'comp'.", call. = FALSE)
      }
    }
    if (length(comp) > 1) {
      if (length(comp) > ncomp) 
        stop("the length of 'comp' must be smaller than or equal to ", 
             ncomp, ".", call. = FALSE)
      if (!is.numeric(comp) || any(comp < 1)) 
        stop("invalid vector for 'comp'.", call. = FALSE)
      if (any(comp > ncomp)) 
        stop("the elements of 'comp' must be smaller or equal than ", 
             ncomp, ".", call. = FALSE)
    }
    comp = round(comp)
    row.names.plot = TRUE
    row.names = mat$names$colnames$X
    if (any(class.object %in% object.pls)) {
      if (all(class(mat) %in% "pls")) {
        keep.X = rep(TRUE, p)
        keep.Y = rep(TRUE, q)
      }
      else {
        keep.X = apply(abs(mat$loadings$X[, comp, drop = FALSE]), 
                       1, sum) > 0
        keep.Y = apply(abs(mat$loadings$Y[, comp, drop = FALSE]), 
                       1, sum) > 0
        row.names = row.names[keep.X]
        col.names = col.names[keep.Y]
      }
      if (mat$mode == "canonical") {
        cord.X = cor(mat$X[, keep.X], mat$variates$X[, 
                                                     comp], use = "pairwise")
        cord.Y = cor(mat$Y[, keep.Y], mat$variates$Y[, 
                                                     comp], use = "pairwise")
      }
      else {
        cord.X = cor(mat$X[, keep.X], mat$variates$X[, 
                                                     comp], use = "pairwise")
        cord.Y = cor(mat$Y[, keep.Y], mat$variates$X[, 
                                                     comp], use = "pairwise")
      }
      mat = cord.X %*% t(cord.Y)
    }
  }
}

smbn<-function(numeric.matrix,netcutoff = NULL){
  Net.edge<-c()
  i=1
  while (i <= ncol(numeric.matrix)){
    tryCatch({
      Xn<-numeric.matrix[,-i]
      Yn<-as.matrix(numeric.matrix[,i])
      colnames(Yn)<-colnames(numeric.matrix[(i)])
      comps<-nrow(numeric.matrix)-(ceiling(nrow(numeric.matrix)/10)+1) #I have an explanation for this line...
      xy.spls <- spls(Xn, Yn, ncomp = comps, mode = "regression")
      
      set.seed(1511)
      nfolds<-nrow(numeric.matrix)
      A<-perf(xy.spls, validation = "Mfold", folds = nfolds, progressBar = FALSE, nrepeat = 100)
      if(all(A$Q2.total[,1]>=0.0975)){
        aux_comp<-1:comps
      }
      else if(all(A$Q2.total[,1]<0.0975)){
        next
      }
      else{
        #j=1
        #while(A$Q2.total[j,1]>=0.0975){ #rule of thumbs of mixOmics...
        #  j=j+1
        #}
        #aux_comp<-1:(j-1)
        aux_comp<-1:min(which(A$Q2.total>=0.0975)) #max or min??????
        
      }
      #print(paste("The used number of components is",length(aux_comp)))
      xy.spls <- spls(Xn, Yn, ncomp = length(aux_comp), mode = "regression")
      net.result<- new_network2(mat=xy.spls, comp = aux_comp, cutoff = netcutoff)  
      
      mat1<-as.matrix(net.result)
      mat2<-as.matrix(mat1[abs(mat1)>=netcutoff,])
      if(nrow(mat2)!=0){
        colnames(mat2)<-colnames(Yn)
        to<-rep(colnames(mat2),dim(mat2)[1])
        mat3<-cbind(to,mat2)
        mat4<-cbind(rownames(mat3),mat3)
        row.names(mat4)<-c(1:(dim(mat3)[1]))
        colnames(mat4)<-c("Source","target","weight")
        Edge.matrix<-mat4
        Net.edge<-rbind(Net.edge,Edge.matrix)
      }
      else{
        stop("There are no edges above network cutoff")
      }
      print(i)
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    i=i+1
  }
  row.names(Net.edge)<-c(1:(dim(Net.edge)[1]))
  return(Net.edge)
}