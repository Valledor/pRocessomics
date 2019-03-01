\name{univariate}
\alias{univariate}
\title{Univariate statistical analysis}

\description{This function performs univariate statistical analyses such as parametrics as t test or ANOVA and non-parametric tests as Wilcox test and Kruskal Wallis, there is also the posibility to adjust obtained pvalues, FDR and post hoc analyses }

\usage{
univariate(datalist,
          initialrow=1,
          initialcolumn=3,
          treatment1col=1,
          treatment2col=2,
          treatment=1,
          parametric=TRUE,
          posthoc=TRUE,
          FDR=TRUE,
          round=5)
}
\arguments{
\item{datalist}{List with different preprocessed omic levels. POL class object.}
\item{initialrow }{First row of numeric data within datasets.}
\item{initialcolumn}{First column of numeric data within datasets.}
\item{treatment1col}{Column number in which treatment 1 is indicated.}
\item{treatment2col}{Column number in which treatment 2 is indicated (if applicable). If only one treatment is set indicate the same column number than treatment1col}
\item{treatment}{1 dataset will be split by treatment 1; 2, dataset will be split by treatment 2; 3, dataset will be split by treatments1&2}
\item{parametric}{logical, this argument indicates whether the employed test should be or not parametric, depending on the number of available treatments if parametric a t test or ANOVA test will be employed, for non parametric data, Wilcox test or kruskal Wallis test will be employed}
\item{posthoc}{logical, this argument indicates whether to do a posthoc test, if parametric=TRUE a TukeyHSD test will be done, for non-parametric analysis Dunn test}
\item{FDR}{logical, if TRUE adjusted pvalues will be provided according to Bonferroni method }
\item{round}{numeric, number of desired decimal digits in the output }
}
\details{
The objective of \code{univariate} function is to perform statistical analyses in order to check variability across treatments . This function requires a \code{POL} object, previously defined by using \code{preprocess_omic_list} and, if desired, \code{transform_varfilter_list}. Two different approaches are implemented within this function, the first related to parametriciy, and the second taking into account the number of treatments to compare. 
}

\value{
List containing mean, standard deviation and selected stats from each omic layer, UNIVAR class object
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Luis Valledor, Laura Lamelas
}

\seealso{\code{\link{transform_varfilter_list}}, \code{\link{feature_selection_list}}, and http://www.miweb.es for more details.}

\examples{
data(xxx)
univariate(datalist = datalist,initialrow = 1,initialcolumn = 5,treatment1col = 1,treatment2col = 1,treatment = 1,omiclevel = NULL,parametric = TRUE,posthoc = FALSE,FDR = TRUE,round = 5)

univariate(datalist = datalist,initialrow = 1,initialcolumn = 5,treatment1col = 1,treatment2col = 1,treatment = 1,omiclevel = NULL,parametric = FALSE,posthoc = TRUE,FDR = TRUE,round = 5)

univariate(datalist = datalist,initialrow = 1,initialcolumn = 5,treatment1col = 1,treatment2col = 1,treatment = 1,omiclevel = NULL,parametric = TRUE,posthoc = TRUE,FDR = TRUE,round = 5)

  }

