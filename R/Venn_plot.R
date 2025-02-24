# Copyright 2019 Luis Valledor / Laura Lamelas
# Last revision 01.2023
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


#' @importFrom grDevices dev.cur dev.off 
#' @importFrom VennDiagram draw.pairwise.venn draw.triple.venn draw.quad.venn 
#' @importFrom graphics pie legend
#' @importFrom plotrix floating.pie 
#' @importFrom gridExtra grid.arrange
#' @importFrom grid gTree textGrob gpar
#' @export
#' @name Venn_plot
#' @title Venn's Diagram plot
#' @description A function to plot the result of a Venn analysis
#' @usage Venn_plot(datalist, Euler.dist = FALSE, alpha=0.5, num.cex=3,font.cex=3)
#' @param datalist An object of class "vennanalysis" provided by Venn_analysis() function
#' @param Euler.dist Boolean indicating whether to draw Euler diagrams when conditions are met or not (Venn Diagrams with movable circles)
#' @param alpha Number giving the alpha transparency of circle areas
#' @param num.cex Number indicating the size of the numbers ploted
#' @param font.cex Number indicating the size of the font depicting the treatments/categories
# @param title Character, plot title
#' @return Venn Diagram plot
#' @author Luis Valledor and Laura Lamelas
#' @details This function will automatically learn how many treatments or categories are to be depicted, in case of treatments > 4 a multi-layer circle plot will be provided


Venn_plot <- function(datalist, Euler.dist = FALSE, alpha = 0.5, num.cex = 3, font.cex = 3) {
  if (!inherits(datalist, "vennanalysis")) stop("A Venn Analysis object is required. Please run Venn_group step first")
  treatments <- datalist[[1]]
  Vennfinalobject <- datalist[[2]]
  
  #Number of treatments
  Treat.No <- length(unique(treatments))
  
  # Custom aesthetics ----
  custom.aes <- Filter(function(x) "color.pRo" %in% class(get(x)), ls(envir= .GlobalEnv))
  if(length(custom.aes) != 0){
    for (i in 1:length(custom.aes)){
      new.aes <- get(custom.aes[i],envir=.GlobalEnv)
      
      switch(class(new.aes)[2], 
             "treatment" = paleta_tratamientos <- new.aes,
             "levels" = paleta_niveles <- new.aes, 
             "continuos" = paleta_continuo <- new.aes,
             "mapman" = paleta_mapman <- new.aes,
             "symbols" = simbolos <- new.aes, 
             "fontsizes" = fontsizes <- new.aes)
    }
    
  }
  
  
  #### Two treatments ####    
  if (Treat.No == 2) {
    a <- as.vector.factor(unique(treatments))
    if (length(which(Vennfinalobject[, 1] == paste(a[1], a[2], sep = "_||_"))) ==
        0) {
      p12 <- 0
    }
    else{
      p12 <-
        Vennfinalobject[which(Vennfinalobject[, 1] == paste(a[1], a[2], sep = "_||_")), 2]
    }
    
    #single treatment groups
    if (length(which(Vennfinalobject[, 1] == a[1])) == 0) {
      area1 <- 0 + p12
    }
    else{
      area1 <- Vennfinalobject[which(Vennfinalobject[, 1] == a[1]), 2] + p12
    }
    
    if (length(which(Vennfinalobject[, 1] == a[2])) == 0) {
      area2 <- 0 + p12
    }
    else{
      area2 <- Vennfinalobject[which(Vennfinalobject[, 1] == a[2]), 2] + p12
    }
    
    
    if (grDevices::dev.cur() > 1)
      grDevices::dev.off()
    myplot <- VennDiagram::draw.pairwise.venn(area1,area2, p12, category = c(a[1], a[2]), euler.d = Euler.dist, scaled = Euler.dist, inverted = FALSE, ext.text = TRUE, ext.percent = rep(0.05, 3), lwd = rep(2, 2), lty = rep("blank", 2), col = rep("black", 2), fill = c("blue", "green"), alpha = rep(alpha, 2), label.col = rep("black", 3), cex= rep(num.cex, 3), fontface = rep("plain", 3), fontfamily = rep("sans", 3), cat.pos = c(0, 0), cat.dist = c(0.05, 0.05), cat.col = c("blue", "green"), cat.cex = rep(font.cex, 2), cat.fontface = rep("plain", 2), cat.fontfamily = rep("sans", 2), cat.just = rep(list(c(0.5, 0.5)), 2), cat.default.pos= "outer", cat.prompts = FALSE, rotation.degree = 0, rotation.centre = c(0.5, 0.5), ind = TRUE, sep.dist = 0.05, offset = 0, cex.prop = NULL, print.mode = "raw", sigdigs = 1)
    
    gridExtra::grid.arrange(grid::gTree(children=myplot), top = grid::textGrob("Two group Venn Diagram", gp=grid::gpar(fontsize=20, fontfamily = "sans", fontfasce="bold")))
    
    class(myplot) <- "Venngrid"
    return(myplot)
  }
  #### Three treatments ####    
  else if (Treat.No == 3) {
    a <- as.vector(unique(treatments))
    
    #triple treatment groups
    if (length(which(Vennfinalobject[, 1] == paste(a[1], a[2], a[3], sep =
                                                   "_||_"))) == 0) {
      n123 <- 0
    }
    else{
      n123 <-
        Vennfinalobject[which(Vennfinalobject[, 1] == paste(a[1], a[2], a[3], sep =
                                                              "_||_")), 2]
    }
    
    #double treatment groups
    if (length(which(Vennfinalobject[, 1] == paste(a[1], a[2], sep = "_||_"))) ==
        0) {
      n12 <- 0 + n123
    }
    else{
      n12 <-
        Vennfinalobject[which(Vennfinalobject[, 1] == paste(a[1], a[2], sep = "_||_")), 2] +
        n123
    }
    
    if (length(which(Vennfinalobject[, 1] == paste(a[2], a[3], sep = "_||_"))) ==
        0) {
      n23 <- 0 + n123
    }
    else{
      n23 <-
        Vennfinalobject[which(Vennfinalobject[, 1] == paste(a[2], a[3], sep = "_||_")), 2] +
        n123
    }
    
    if (length(which(Vennfinalobject[, 1] == paste(a[1], a[3], sep = "_||_"))) ==
        0) {
      n13 <- 0 + n123
    }
    else{
      n13 <-
        Vennfinalobject[which(Vennfinalobject[, 1] == paste(a[1], a[3], sep = "_||_")), 2] +
        n123
    }
    
    #single treatment groups
    if (length(which(Vennfinalobject[, 1] == a[1])) == 0) {
      area1 <- 0 + n12 + n13 - n123
    }
    else{
      area1 <-
        Vennfinalobject[which(Vennfinalobject[, 1] == a[1]), 2] + n12 + n13 - n123
    }
    
    if (length(which(Vennfinalobject[, 1] == a[2])) == 0) {
      area2 <- 0 + n12 + n23 - n123
    }
    else{
      area2 <-
        Vennfinalobject[which(Vennfinalobject[, 1] == a[2]), 2] + n12 + n23 - n123
    }
    
    if (length(which(Vennfinalobject[, 1] == a[3])) == 0) {
      area3 <- 0 + n13 + n23 - n123
    }
    else{
      area3 <-
        Vennfinalobject[which(Vennfinalobject[, 1] == a[3]), 2] + n13 + n23 - n123
    }
    
    #he aqui el plot... basico
    if (grDevices::dev.cur() > 1){
      grDevices::dev.off()}
    myplot <- VennDiagram::draw.triple.venn(area1,area2, area3,n12, n23, n13, n123, category =c(a[1], a[2], a[3]),
                                            rotation = 1,reverse = FALSE, euler.d = Euler.dist, scaled = Euler.dist, lwd = rep(2, 3),
                                            lty = rep("blank", 3), col = rep("black", 3), fill = paleta_tratamientos[1:3], alpha = rep(alpha, 3),
                                            label.col = rep("black", 7), cex= rep(num.cex, 7), fontface = rep("plain", 7), fontfamily = rep("sans", 7),
                                            cat.pos = c(-40, 40, 180),cat.dist = c(0.05, 0.05, 0.05), cat.col = paleta_tratamientos[1:3],
                                            cat.cex = rep(font.cex, 3), cat.fontface = rep("plain", 3), cat.fontfamily = rep("sans", 3),
                                            cat.just = list(c(0.5, 1), c(0.5, 1), c(0.5, 0.25)), cat.default.pos = "outer", cat.prompts = FALSE,
                                            rotation.degree = 0, rotation.centre = c(0.5, 0.5),ind = TRUE,sep.dist = 0.05,offset = 0,cex.prop = NULL,
                                            print.mode = "raw", sigdigs = 1)
    gridExtra::grid.arrange(grid::gTree(children=myplot), top = grid::textGrob("Three group Venn Diagram", gp=grid::gpar(fontsize=20, fontfamily = "sans", fontfasce="bold")))
    class(myplot) <- "Venngrid"
    return(myplot)
  }
  #### Four treatments ####
  else if (Treat.No == 4) {
    #require(VennDiagram)
    aa <- as.vector(unique(treatments))
    #quadruple treatment groups
    if (length(which(Vennfinalobject[, 1] == paste(aa[1], aa[2], aa[3], aa[4], sep =
                                                   "_||_"))) == 0) {
      m1234 <- 0
    }
    else{
      m1234 <-
        Vennfinalobject[which(Vennfinalobject[, 1] == paste(aa[1], aa[2], aa[3], aa[4], sep =
                                                              "_||_")), 2]
    }
    
    #triple treatment groups
    if (length(which(Vennfinalobject[, 1] == paste(aa[1], aa[2], aa[3], sep =
                                                   "_||_"))) == 0) {
      m123 <- 0 + m1234
    }
    else{
      m123 <-
        Vennfinalobject[which(Vennfinalobject[, 1] == paste(aa[1], aa[2], aa[3], sep =
                                                              "_||_")), 2] + m1234
    }
    
    if (length(which(Vennfinalobject[, 1] == paste(aa[1], aa[3], aa[4], sep =
                                                   "_||_"))) == 0) {
      m134 <- 0 + m1234
    }
    else{
      m134 <-
        Vennfinalobject[which(Vennfinalobject[, 1] == paste(aa[1], aa[3], aa[4], sep =
                                                              "_||_")), 2] + m1234
    }
    
    if (length(which(Vennfinalobject[, 1] == paste(aa[1], aa[2], aa[4], sep =
                                                   "_||_"))) == 0) {
      m124 <- 0 + m1234
    }
    else{
      m124 <-
        Vennfinalobject[which(Vennfinalobject[, 1] == paste(aa[1], aa[2], aa[4], sep =
                                                              "_||_")), 2] + m1234
    }
    
    if (length(which(Vennfinalobject[, 1] == paste(aa[2], aa[3], aa[4], sep = "_||_"))) ==
        0) {
      m234 <- 0 + m1234
    }
    else{
      m234 <-
        Vennfinalobject[which(Vennfinalobject[, 1] == paste(aa[2], aa[3], aa[4], sep =
                                                              "_||_")), 2] + m1234
    }
    
    #double treatment groups
    if (length(which(Vennfinalobject[, 1] == paste(aa[1], aa[2], sep = "_||_"))) ==
        0) {
      m12 <- 0 + m123 + m124 - m1234
    }
    else{
      m12 <-
        Vennfinalobject[which(Vennfinalobject[, 1] == paste(aa[1], aa[2], sep =
                                                              "_||_")), 2] + m123 + m124 - m1234
    }
    
    if (length(which(Vennfinalobject[, 1] == paste(aa[1], aa[3], sep = "_||_"))) ==
        0) {
      m13 <- 0 + m123 + m134 - m1234
    }
    else{
      m13 <-
        Vennfinalobject[which(Vennfinalobject[, 1] == paste(aa[1], aa[3], sep =
                                                              "_||_")), 2] + m123 + m134 - m1234
    }
    
    if (length(which(Vennfinalobject[, 1] == paste(aa[1], aa[4], sep = "_||_"))) ==
        0) {
      m14 <- 0 + m124 + m134 - m1234
    }
    else{
      m14 <-
        Vennfinalobject[which(Vennfinalobject[, 1] == paste(aa[1], aa[4], sep =
                                                              "_||_")), 2] + m124 + m134 - m1234
    }
    
    if (length(which(Vennfinalobject[, 1] == paste(aa[2], aa[3], sep = "_||_"))) ==
        0) {
      m23 <- 0 + m123 + m234 - m1234
    }
    else{
      m23 <-
        Vennfinalobject[which(Vennfinalobject[, 1] == paste(aa[2], aa[3], sep =
                                                              "_||_")), 2] + m123 + m234 - m1234
    }
    
    if (length(which(Vennfinalobject[, 1] == paste(aa[2], aa[4], sep = "_||_"))) ==
        0) {
      m24 <- 0 + m124 + m234 - m1234
    }
    else{
      m24 <-
        Vennfinalobject[which(Vennfinalobject[, 1] == paste(aa[2], aa[4], sep =
                                                              "_||_")), 2] + m124 + m234 - m1234
    }
    
    if (length(which(Vennfinalobject[, 1] == paste(aa[3], aa[4], sep = "_||_"))) ==
        0) {
      m34 <- 0 + m134 + m234 - m1234
    }
    else{
      m34 <-
        Vennfinalobject[which(Vennfinalobject[, 1] == paste(aa[3], aa[4], sep =
                                                              "_||_")), 2] + m134 + m234 - m1234
    }
    
    
    #single treatment groups
    if (length(which(Vennfinalobject[, 1] == aa[1])) == 0) {
      aream1 <- 0 + m12 + m13 + m14 - m123 - m124 - m134 + m1234
    }
    else{
      aream1 <-
        Vennfinalobject[which(Vennfinalobject[, 1] == aa[1]), 2] + m12 + m13 + m14 -
        m123 - m124 - m134 + m1234
    }
    
    if (length(which(Vennfinalobject[, 1] == aa[2])) == 0) {
      aream2 <- 0 + m12 + m23 + m24 - m123 - m124 - m234 + m1234
    }
    else{
      aream2 <-
        Vennfinalobject[which(Vennfinalobject[, 1] == aa[2]), 2] + m12 + m23 + m24 -
        m123 - m124 - m234 + m1234
    }
    
    if (length(which(Vennfinalobject[, 1] == aa[3])) == 0) {
      aream3 <- 0 + m13 + m23 + m34 - m123 - m134 - m234 + m1234
    }
    else{
      aream3 <-
        Vennfinalobject[which(Vennfinalobject[, 1] == aa[3]), 2] + m13 + m23 + m34 -
        m123 - m134 - m234 + m1234
    }
    
    if (length(which(Vennfinalobject[, 1] == aa[4])) == 0) {
      aream4 <- 0 + m14 + m24 + m34 - m124 - m134 - m234 + m1234
    }
    else{
      aream4 <-
        Vennfinalobject[which(Vennfinalobject[, 1] == aa[4]), 2] + m14 + m24 + m34 -
        m124 - m134 - m234 + m1234
    }
    
    if (grDevices::dev.cur() > 1)
      grDevices::dev.off()
    myplot <- VennDiagram::draw.quad.venn(aream1, aream2, aream3, aream4, m12, m13, m14, m23, m24, m34, m123, m124, m134, m234, m1234, category = aa, lwd = rep(2, 4), lty = rep("blank", 4), 
                                          col = rep("black", 4), fill = paleta_tratamientos[1:4], alpha = rep(alpha, 4), label.col = rep("black", 15), cex = rep(num.cex, 15), 
                                          fontface = rep("plain", 15), fontfamily = rep("sans", 15), cat.pos = c(-15, 15, 0, 0), cat.dist = c(0.22, 0.22, 0.11, 0.11), 
                                          cat.col = paleta_tratamientos[1:4], cat.cex = rep(font.cex, 4), cat.fontface = rep("plain", 4), cat.fontfamily = rep("sans", 4), 
                                          cat.just = rep(list(c(0.5, 0.5)), 4), rotation.degree = 0, rotation.centre = c(0.5, 0.5), ind = TRUE, cex.prop = NULL, print.mode = "raw", sigdigs = 3, 
                                          direct.area = FALSE, area.vector = 0)
    
    gridExtra::grid.arrange(grid::gTree(children=myplot), top = grid::textGrob("Four group Venn Diagram", gp=grid::gpar(fontsize=20, fontfamily = "sans", fontfasce="bold")))
    class(myplot) <- "Venngrid"
    return(myplot)
  }
  #### Sunburst ####
  else if (Treat.No >= 5) {
    vennmat <- datalist[[2]]
    treatments <- datalist[[1]]
    sections <- dim(vennmat)[1]
    #ahora genero un objeto, que sera el tamano,es decir, la frecuencia de cada section
    sizesection <- vennmat[, 2]
    #lo siguiente es generar el objeto que me diga que tratamiento esta en cada sector, un boolean T/F, que luego se convirtio en una matriz de 0 y 1
    #va a tener que ser un bucle
    Treats <- as.vector.factor(unique(treatments))
    present.object <- c()
    
    for (k in 1:(length(Treats))) {
      is.there <- c()
      is.present <- c()
      for (i in 1:sections) {
        is.there <- grepl(Treats[k], ((datalist[[2]])[i, 1]))
        if (is.there == T) {
          is.there <- 1
        } else{
          is.there <- 0
        }
        is.there <- is.there * sizesection[i]
        is.present <- rbind(is.present, is.there)
      }
      present.object <- cbind(present.object, is.present)
    }
    colnames(present.object) <- Treats
    rownames(present.object) <- vennmat[, 1]
    present.object <- as.data.frame(present.object)
    coloursvec <- paleta_tratamientos[1:Treat.No]
    
    colourmat <- c()
    for (k in 1:(length(Treats))) {
      is.on <- c()
      is.off <- c()
      for (i in 1:sections) {
        if (present.object[i, k] == 0) {
          is.on <- "white"
        } else{
          is.on <- coloursvec[k]
        }
        is.off <- rbind(is.off, is.on)
      }
      colourmat <- cbind(colourmat, is.off)
    }
    
    colourmat <- as.data.frame(colourmat)
    colnames(colourmat) <- Treats
    rownames(colourmat) <- vennmat[, 1]
    iniR <- 0.15
    xsect <- as.vector(sizesection)
    graphics::par(mai=c(1,1,1,1))
    graphics::pie(0, 0, x = xsect, edges = 1000, radius = (Treat.No + 1) * iniR, col = as.vector(colourmat[, Treat.No]), border = NA, labels = '')
    
    for (i in ((length(Treats)) - 1):1) {
     plotrix::floating.pie(0, 0, x = xsect, edges = 1000, radius = (i + 1) * iniR, col = as.vector(colourmat[, i]), border = NA)
    }
    
    plotrix::floating.pie(0, 0, x = 1, edges = 1000, radius = iniR, col = "white", border = NA)
    graphics::legend(x = -1.5/2,y = 1.1, Treats,col = as.character(coloursvec),pch = 19, bty = 'n',ncol = length(Treats))
    graphics::title("Sunburst plot")
    #graphics::title(paste("Sunburst plot of ",deparse(substitute(datalist))))
    
   invisible(a<-list("sunburst",datalist))
   
  }
  
}
