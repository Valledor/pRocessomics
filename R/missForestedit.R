#
# Copyright 2019 Luis Valledor / Laura Lamelas
# Last revision 04.11.2019
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

#' @importFrom foreach getDoParWorkers foreach %dopar%
#' @importFrom itertools isplitVector
#' @importFrom randomForest randomForest
#' @importFrom stats predict var
#' @importFrom iterators idiv
#' @importFrom missForest mixError 

missForestedit<-function (xmis, maxiter = 10, ntree = 100, variablewise = FALSE, decreasing = FALSE, verbose = FALSE, mtry = floor(sqrt(ncol(xmis))), replace = TRUE, classwt = NULL, cutoff = NULL, strata = NULL, 
                          sampsize = NULL, nodesize = NULL, maxnodes = NULL, xtrue = NA, parallelize = c("no", "variables", "forests")) {
  suppressPackageStartupMessages(requireNamespace("randomForest"))
  n <- nrow(xmis)
  p <- ncol(xmis)
  if (!is.null(classwt)) 
    stopifnot(length(classwt) == p, typeof(classwt) == "list")
  if (!is.null(cutoff)) 
    stopifnot(length(cutoff) == p, typeof(cutoff) == "list")
  if (!is.null(strata)) 
    stopifnot(length(strata) == p, typeof(strata) == "list")
  if (!is.null(nodesize)) 
    stopifnot(length(nodesize) == 2)
  if (any(apply(is.na(xmis), 2, sum) == n)) {
    indCmis <- which(apply(is.na(xmis), 2, sum) == n)
    xmis <- xmis[, -indCmis]
    p <- ncol(xmis)
    cat("  removed variable(s)", indCmis, "due to the missingness of all entries\n")
  }
  parallelize <- match.arg(parallelize)
  if (parallelize %in% c("variables", "forests")) {
    if (foreach::getDoParWorkers() == 1) {
      stop("You must register a 'foreach' parallel backend to run 'missForest' in parallel. Set 'parallelize' to 'no' to compute serially.")
    }
    else if (verbose) {
      if (parallelize == "variables") {
        cat("  parallelizing over the variables of the input data matrix 'xmis'\n")
      }
      else {
        cat("  parallelizing computation of the random forest model objects\n")
      }
    }
    if (foreach::getDoParWorkers() > p) {
      stop("The number of parallel cores should not exceed the number of variables (p=", 
           p, ")")
    }
  }
  ximp <- xmis
  xAttrib <- lapply(xmis, attributes)
  varType <- character(p)
  for (t.co in 1:p) {
    if (is.null(xAttrib[[t.co]])) {
      varType[t.co] <- "numeric"
      ximp[is.na(xmis[, t.co]), t.co] <- mean(xmis[, t.co], 
                                              na.rm = TRUE)
    }
    else {
      varType[t.co] <- "factor"
      max.level <- max(table(ximp[, t.co]))
      class.assign <- sample(names(which(max.level == summary(ximp[, 
                                                                   t.co]))), 1)
      if (class.assign != "NA's") {
        ximp[is.na(xmis[, t.co]), t.co] <- class.assign
      }
      else {
        while (class.assign == "NA's") {
          class.assign <- sample(names(which(max.level == 
                                               summary(ximp[, t.co]))), 1)
        }
        ximp[is.na(xmis[, t.co]), t.co] <- class.assign
      }
    }
  }
  NAloc <- is.na(xmis)
  noNAvar <- apply(NAloc, 2, sum)
  sort.j <- order(noNAvar)
  if (decreasing) 
    sort.j <- rev(sort.j)
  sort.noNAvar <- noNAvar[sort.j]
  nzsort.j <- sort.j[sort.noNAvar > 0]
  if (parallelize == "variables") {
    "%cols%" <- get("%dopar%")
    idxList <- as.list(itertools::isplitVector(nzsort.j, chunkSize = foreach::getDoParWorkers()))
  }
  Ximp <- vector("list", maxiter)
  iter <- 0
  k <- length(unique(varType))
  convNew <- rep(0, k)
  convOld <- rep(Inf, k)
  OOBerror <- numeric(p)
  names(OOBerror) <- varType
  if (k == 1) {
    if (unique(varType) == "numeric") {
      names(convNew) <- c("numeric")
    }
    else {
      names(convNew) <- c("factor")
    }
    convergence <- c()
    OOBerr <- numeric(1)
  }
  else {
    names(convNew) <- c("numeric", "factor")
    convergence <- matrix(NA, ncol = 2)
    OOBerr <- numeric(2)
  }
  stopCriterion <- function(varType, convNew, convOld, iter, 
                            maxiter) {
    k <- length(unique(varType))
    if (k == 1) {
      (convNew < convOld) & (iter < maxiter)
    }
    else {
      ((convNew[1] < convOld[1]) | (convNew[2] < convOld[2])) & 
        (iter < maxiter)
    }
  }
  while (stopCriterion(varType, convNew, convOld, iter, maxiter)) {
    if (iter != 0) {
      convOld <- convNew
      OOBerrOld <- OOBerr
    }
    cat("  missForest iteration", iter + 1, "in progress...")
    t.start <- proc.time()
    ximp.old <- ximp
    if (parallelize == "variables") {
      for (idx in idxList) {
        
        results <- foreach::foreach(varInd = idx, .packages = "randomForest") %cols% 
        {
          obsi <- !NAloc[, varInd]
          misi <- NAloc[, varInd]
          obsY <- ximp[obsi, varInd]
          obsX <- ximp[obsi, seq(1, p)[-varInd]]
          misX <- ximp[misi, seq(1, p)[-varInd]]
          typeY <- varType[varInd]
          if (typeY == "numeric") {
            set.seed(1511)
            RF <- randomForest::randomForest(x = obsX, y = obsY, 
                                             ntree = ntree, mtry = mtry, replace = replace, 
                                             sampsize = if (!is.null(sampsize)) 
                                               sampsize[[varInd]]
                                             else if (replace) 
                                               nrow(obsX)
                                             else ceiling(0.632 * nrow(obsX)), nodesize = if (!is.null(nodesize)) 
                                               nodesize[1]
                                             else 1, maxnodes = if (!is.null(maxnodes)) 
                                               maxnodes
                                             else NULL)
            oerr <- RF$mse[ntree]
            misY <- stats::predict(RF, misX)
          }
          else {
            obsY <- factor(obsY)
            summarY <- summary(obsY)
            if (length(summarY) == 1) {
              oerr <- 0
              misY <- factor(rep(names(summarY), length(misi)))
            }
            else {
              set.seed(1511)
              RF <- randomForest::randomForest(x = obsX, y = obsY, 
                                               ntree = ntree, mtry = mtry, replace = replace, 
                                               classwt = if (!is.null(classwt)) 
                                                 classwt[[varInd]]
                                               else rep(1, nlevels(obsY)), cutoff = if (!is.null(cutoff)) 
                                                 cutoff[[varInd]]
                                               else rep(1/nlevels(obsY), nlevels(obsY)), 
                                               strata = if (!is.null(strata)) 
                                                 strata[[varInd]]
                                               else obsY, sampsize = if (!is.null(sampsize)) 
                                                 sampsize[[varInd]]
                                               else if (replace) 
                                                 nrow(obsX)
                                               else ceiling(0.632 * nrow(obsX)), nodesize = if (!is.null(nodesize)) 
                                                 nodesize[2]
                                               else 5, maxnodes = if (!is.null(maxnodes)) 
                                                 maxnodes
                                               else NULL)
              oerr <- RF$err.rate[[ntree, 1]]
              misY <- stats::predict(RF, misX)
            }
          }
          list(varInd = varInd, misY = misY, oerr = oerr)
        }
        for (res in results) {
          misi <- NAloc[, res$varInd]
          ximp[misi, res$varInd] <- res$misY
          OOBerror[res$varInd] <- res$oerr
        }
      }
    }
    else {
      for (s in 1:p) {
        varInd <- sort.j[s]
        if (noNAvar[[varInd]] != 0) {
          obsi <- !NAloc[, varInd]
          misi <- NAloc[, varInd]
          obsY <- ximp[obsi, varInd]
          obsX <- ximp[obsi, seq(1, p)[-varInd]]
          misX <- ximp[misi, seq(1, p)[-varInd]]
          typeY <- varType[varInd]
          if (typeY == "numeric") {
            if (parallelize == "forests") {
              xntree <- NULL
              set.seed(1511)
              RF <- foreach::foreach(xntree = iterators::idiv(ntree, chunks = foreach::getDoParWorkers()), 
                                     .combine = "combine", .multicombine = TRUE, 
                                     .packages = "randomForest") %dopar% {
                                       randomForest::randomForest(x = obsX, y = obsY, ntree = xntree, 
                                                                  mtry = mtry, replace = replace, sampsize = if (!is.null(sampsize)) 
                                                                    sampsize[[varInd]]
                                                                  else if (replace) 
                                                                    nrow(obsX)
                                                                  else ceiling(0.632 * nrow(obsX)), nodesize = if (!is.null(nodesize)) 
                                                                    nodesize[1]
                                                                  else 1, maxnodes = if (!is.null(maxnodes)) 
                                                                    maxnodes
                                                                  else NULL)
                                     }
              OOBerror[varInd] <- mean((stats::predict(RF) - 
                                          RF$y)^2, na.rm = TRUE)
            }
            else {
              set.seed(1511)        
              RF <- randomForest::randomForest(x = obsX, y = obsY, 
                                               ntree = ntree, mtry = mtry, replace = replace, 
                                               sampsize = if (!is.null(sampsize)) 
                                                 sampsize[[varInd]]
                                               else if (replace) 
                                                 nrow(obsX)
                                               else ceiling(0.632 * nrow(obsX)), nodesize = if (!is.null(nodesize)) 
                                                 nodesize[1]
                                               else 1, maxnodes = if (!is.null(maxnodes)) 
                                                 maxnodes
                                               else NULL)
              OOBerror[varInd] <- RF$mse[ntree]
            }
            misY <- stats::predict(RF, misX)
          }
          else {
            obsY <- factor(obsY)
            summarY <- summary(obsY)
            if (length(summarY) == 1) {
              misY <- factor(rep(names(summarY), sum(misi)))
            }
            else {
              if (parallelize == "forests") {
                set.seed(1511)
                RF <- foreach::foreach(xntree = iterators::idiv(ntree, chunks = foreach::getDoParWorkers()), 
                                       .combine = "combine", .multicombine = TRUE, 
                                       .packages = "randomForest") %dopar% 
                                       {
                                         set.seed(1511)
                                         randomForest::randomForest(x = obsX, y = obsY, 
                                                                    ntree = xntree, mtry = mtry, replace = replace, 
                                                                    classwt = if (!is.null(classwt)) 
                                                                      classwt[[varInd]]
                                                                    else rep(1, nlevels(obsY)), cutoff = if (!is.null(cutoff)) 
                                                                      cutoff[[varInd]]
                                                                    else rep(1/nlevels(obsY), nlevels(obsY)), 
                                                                    strata = if (!is.null(strata)) 
                                                                      strata[[varInd]]
                                                                    else obsY, sampsize = if (!is.null(sampsize)) 
                                                                      sampsize[[varInd]]
                                                                    else if (replace) 
                                                                      nrow(obsX)
                                                                    else ceiling(0.632 * nrow(obsX)), 
                                                                    nodesize = if (!is.null(nodesize)) 
                                                                      nodesize[2]
                                                                    else 5, maxnodes = if (!is.null(maxnodes)) 
                                                                      maxnodes
                                                                    else NULL)
                                       }
                ne <- as.integer(stats::predict(RF)) != as.integer(RF$y)
                ne <- ne[!is.na(ne)]
                OOBerror[varInd] <- sum(ne)/length(ne)
              }
              else {
                set.seed(1511)
                RF <- randomForest::randomForest(x = obsX, y = obsY, 
                                                 ntree = ntree, mtry = mtry, replace = replace, 
                                                 classwt = if (!is.null(classwt)) 
                                                   classwt[[varInd]]
                                                 else rep(1, nlevels(obsY)), cutoff = if (!is.null(cutoff)) 
                                                   cutoff[[varInd]]
                                                 else rep(1/nlevels(obsY), nlevels(obsY)), 
                                                 strata = if (!is.null(strata)) 
                                                   strata[[varInd]]
                                                 else obsY, sampsize = if (!is.null(sampsize)) 
                                                   sampsize[[varInd]]
                                                 else if (replace) 
                                                   nrow(obsX)
                                                 else ceiling(0.632 * nrow(obsX)), nodesize = if (!is.null(nodesize)) 
                                                   nodesize[2]
                                                 else 5, maxnodes = if (!is.null(maxnodes)) 
                                                   maxnodes
                                                 else NULL)
                OOBerror[varInd] <- RF$err.rate[[ntree, 
                                                 1]]
              }
              misY <- stats::predict(RF, misX)
            }
          }
          ximp[misi, varInd] <- misY
        }
      }
    }
    cat("done!\n")
    iter <- iter + 1
    Ximp[[iter]] <- ximp
    t.co2 <- 1
    for (t.type in names(convNew)) {
      t.ind <- which(varType == t.type)
      if (t.type == "numeric") {
        convNew[t.co2] <- sum((ximp[, t.ind] - ximp.old[, 
                                                        t.ind])^2)/sum(ximp[, t.ind]^2)
      }
      else {
        dist <- sum(as.character(as.matrix(ximp[, t.ind])) != 
                      as.character(as.matrix(ximp.old[, t.ind])))
        convNew[t.co2] <- dist/(n * sum(varType == "factor"))
      }
      t.co2 <- t.co2 + 1
    }
    if (!variablewise) {
      NRMSE <- sqrt(mean(OOBerror[varType == "numeric"])/stats::var(as.vector(as.matrix(xmis[, 
                                                                                             varType == "numeric"])), na.rm = TRUE))
      PFC <- mean(OOBerror[varType == "factor"])
      if (k == 1) {
        if (unique(varType) == "numeric") {
          OOBerr <- NRMSE
          names(OOBerr) <- "NRMSE"
        }
        else {
          OOBerr <- PFC
          names(OOBerr) <- "PFC"
        }
      }
      else {
        OOBerr <- c(NRMSE, PFC)
        names(OOBerr) <- c("NRMSE", "PFC")
      }
    }
    else {
      OOBerr <- OOBerror
      names(OOBerr)[varType == "numeric"] <- "MSE"
      names(OOBerr)[varType == "factor"] <- "PFC"
    }
    if (any(!is.na(xtrue))) {
      err <- suppressWarnings(missForest::mixError(ximp, xmis, xtrue))
    }
    if (verbose) {
      delta.start <- proc.time() - t.start
      if (any(!is.na(xtrue))) {
        cat("    error(s):", err, "\n")
      }
      cat("    estimated error(s):", OOBerr, "\n")
      cat("    difference(s):", convNew, "\n")
      cat("    time:", delta.start[3], "seconds\n\n")
    }
  }
  if (iter == maxiter) {
    if (any(is.na(xtrue))) {
      out <- list(ximp = Ximp[[iter]], OOBerror = OOBerr)
    }
    else {
      out <- list(ximp = Ximp[[iter]], OOBerror = OOBerr, 
                  error = err)
    }
  }
  else {
    if (any(is.na(xtrue))) {
      out <- list(ximp = Ximp[[iter - 1]], OOBerror = OOBerrOld)
    }
    else {
      out <- list(ximp = Ximp[[iter - 1]], OOBerror = OOBerrOld, 
                  error = suppressWarnings(missForest::mixError(Ximp[[iter - 
                                                                        1]], xmis, xtrue)))
    }
  }
  class(out) <- "missForest"
  return(out)
}
