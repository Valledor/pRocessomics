\name{transform_varfilter_list}
\alias{transform_varfilter_list}
\title{Transformation and Filtering of Variables of Omic Datasets}

\description{
This function processes a POL object with the different omic levels to be analized in order to transform abundance values and select variables based on consistency, significant values, or coefficient of variation.
}

\usage{
transform_varfilter_list(datalist,
                    initialrow=1,
                    initialcolumn=3,
                    treatment1col = 1,
                    treatment2col = 2,
                    treatment = 3,
                    transf=NULL,
                    varsel=FALSE,
                    varselthld=0.25,
                    varcoef=FALSE,
                    varcoefthld=75)
                    }
\arguments{
\item{datalist}{A "POL" object.}
\item{initialrow }{First row of numeric data within datasets.}
\item{initialcolumn}{First column of numeric data within datasets.}
\item{treatment1col}{Column number in which treatment 1 is indicated.}
\item{treatment2col}{Column number in which treatment 2 is indicated (if applicable). If only one treatment is set indicate the same column number than treatment1col}
\item{treatment}{1 dataset will be split by treatment 1; 2, dataset will be split by treatment 2; 3, dataset will be split by treatments1&2}
\item{transf}{Method for transforming variables (z, z-centered, Log10, Log10plus1, box-cox-lambda, sqrt, cubic, none; see details) This can be provided as a single value or a vector with the same number of elements than list indicating specific transformation methods for each omic level.
}
\item{varsel}{boleean. If TRUE variables will be filtered by the number of significant values.This can be provided as a single value or a vector with the same number of elements than list to indicating wich omic sets must be processed.}
\item{varselthld}{Indicates the mininum number in \% of significant values required to consider variables. This can be provided as a single value or a vector with the same number of elements than list indicating specific thresholds.}
\item{varcoef}{boleean. If TRUE variables will also be filtered by coefficient of variation.}
\item{varcoefthld}{Indicates the coefficient of variation above which variables will be filtered out. This can be provided as a single value or a vector with the same number of elements than list indicating specific thresholds.}
}

\details{
The objective of \code{transform_varfilter_list} is transforming dataset values and removing those inconsistent or noisy variables. Transformation of values can be performed employing different methods z-score, z, z-score centered, z-centered, logarithm base 10, Log10, logarithm base 10 after adding 1 to all values (avoid inf), Log10plus1, box-cox-lambda transform, box-cox-lambda, square root, sqrt, cubic root, cubic, or none. The selection of consistent variables is performed based on the number of significant values available for each variable (varselect=TRUE) or by its coefficient of variation (varcoef=TRUE).
}

\value{
  Returns a "POL" object, a list with the filtered and/or transformed dataset.
}

\author{
  Luis Valledor, Laura Lamelas
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
