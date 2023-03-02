#' @name univariate_plot
#' @title Univariate statistical analysis plot
#' @description This function plots univariate statistical analyses by showing the p value distribution
#' @usage univariate_plot(univarobject, plottype = "distribution", annotatefile = NULL, 
#' savetofile = F, filename = NULL)
#' @param univarobject "UNIVAR" class object from univariate_analysis function
#' @param plottype 
#' \itemize{
#' \item distribution: scatter plot of p-values
#' \item histogram}
#' @param annotatefile optional, pRoAnnot class object created with pRocessomics::importannotation function
#' @param savetofile TRUE or FALSE, indicating if the plots should be 
#' @param filename only applicable if savetofile set to TRUE

#' @return p-value distribution plot or histogram
#' @author Luis Valledor and Laura Lamelas
#' @seealso univariate, univariate_wizard
#' 
#' @export

univariate_plot<-function(univarobject, plottype = "distribution",annotatefile=NULL,savetofile=F,filename=NULL){
  if(plottype %in% c("distribution","histogram") == F) stop("\nPlease select a valid plottype. Available options are distribution and histogram")
  objname<-deparse(substitute(univarobject))
  
  if(!inherits(univarobject, "UNIVAR")) stop("\nAn UNIVAR object required. Please run univariate function first")
  
  omicstats<-univarobject$composite
  if(!is.null(annotatefile)){
  omicstats2<-lapply(omicstats, function(x) {
    as.data.frame(x)
    x[,c("p-value","DESCRIPTION"),drop=F] 
  })}
  else{
    omicstats2<-lapply(omicstats, function(x) {
      as.data.frame(x)
      x[,c("p-value"),drop=F] 
    })}
  if(plottype=="histogram") plotcases <- TRUE
  if(plottype=="distribution") plotcases <- FALSE
 plots<- lapply(names(omicstats2), function(x) statscatter(omicstats2[[x]], x, savetofile,filename, plotcases))
}




