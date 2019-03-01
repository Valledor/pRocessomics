\name{feature_selection_list}
\alias{feature_selection_list}
\title{Feature Selection in Omic Datasets}

\description{
This function processes a POL object with the different omic levels to be analized in order to select featured variables.
}

\usage{
feature_selection_list(datalist,
                    initialrow=1,
                    initialcolumn=3,
                    treatment1col=1,
                    treatment2col=2,
                    treatment=1,
                    method=NULL,
                    threshold=NULL,
                    parametric=TRUE,
                    stat="p")
                    }
\arguments{
\item{datalist}{A "POL" object.}
\item{initialrow }{First row of numeric data within datasets.}
\item{initialcolumn}{First column of numeric data within datasets.}
\item{treatment1col}{Column number in which treatment 1 is indicated.}
\item{treatment2col}{Column number in which treatment 2 is indicated (if applicable). If only one treatment is set indicate the same column number than treatment1col}
\item{treatment}{1 dataset will be split by treatment 1; 2, dataset will be split by treatment 2; 3, dataset will be split by treatments1&2}
\item{method}{Method for selecting variables: "IQR"- interquartilic range; "Stat" - ANOVA/Kruskal Wallis.}
\item{threshold}{For "IQR": \% over the average IQR; For "Stat": p or q value}
\item{parametric}{boleean. If TRUE ANOVA test will be used, if FALSE Kruskal Wallis.}
\item{stat}{Indicates which statistic will be used p-values,"p", or q-values, "q". }
}

\details{
The objective of \code{feature_selection_list} is removing those variables who do not greatly change between treatments, known to later introduce noise in multivariate analyses reducing the power of these approaches. This function requires a \code{POL} object, previously defined by using \code{preprocess_omic_list} and, if desired, \code{transform_varfilter_list}. Two different approaches are implemented within this function: first one is based on selecting those variables that have intercuartilic ranges, IQR, greater than average within all studied variables. In this case the \code{threshold} indicates the \% above average IQR and only variables with greater or equal IQR will be selected. The second strategy is based on applying ANOVA or Kruskal Wallis and selected those variables exhibiting a p or q value smaller or equal than \code{threshold}.
}

\value{
  Returns a "POL" object, a list with the featured dataset.
}

\author{
  Luis Valledor
}

\references{
  \code{\link{citasd}}.
}


\seealso{
\code{\link{preprocess_omic_list}},
\code{\link{transform_varfilter_list}},
and http://www.miweb.es for more details.
}

\examples{
data(xxx)
my_processed_list <- preprocess_omic_list(my_list, initialrow = 1, initialcolumn = 3, treatment1col = 1, treatment2col = 2, treatment = 3, imputation = "RF", abdbal="AvgIntensity", threshold=0.35, k=3, parallel=TRUE)
}
