\name{da_analysis}
\alias{da_analysis}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{A function to perform a sparce partial least squared discriminant analysis}

\description{This function takes a preprocessed omic datalist and performs a multivariate spls-da analysis}

\usage{
da_analysis(datalist,
                    annotation=FALSE,
                    initialrow=1,
                    initialcolumn=3,
                    treatment1col=1,
                    treatment2col=2,
                    treatment=1,
                    omiclevel="all",
                    KeepX=NULL,
                    autotune=TRUE,
                    performance=FALSE
                    folds=10)
                    }

\arguments{
  \item{datalist}{input preprocessed values, a POL object is required}
  \item{annotation}{data rame containing identifiers, descriptions and bincodes from Mapman annotation}
  \item{initialrow}{First row of numeric data within datasets.}
  \item{initialcolumn}{First column of numeric data within datasets.}
  \item{treatment1col}{Column number in which treatment 1 is indicated.}
  \item{treatment2col}{Column number in which treatment 2 is indicated (if applicable). If only one treatment is set indicate the same column number than treatment1col}
  \item{treatment}{1 dataset will be split by treatment 1; 2, dataset will be split by treatment 2; 3, dataset will be split by treatments1&2}
  \item{omiclevel}{name of the omic layer to analyse}
  \item{KeepX}{Number of variables to take into account for the analysis, if set to NULL all variables will enter in the analysis}
  \item{autotune}{logical, if set to TRUE several KeepX values will be tested and the one with the less error will be choose}
  \item{perfomance}{logical, analysis of the generated model performance when one and two components are used }
  \item{folds}{the folds in the Mfold cross-validation, maximum allowed =10}
}
\details{The aim of this function is to perform a supervised analysis combining the assayed experimental treatments against elemant values such as protein/metabolites/transcrips quantifications}
\value{da_analysis returns an object of class "daanalysis", a list that contains the following components:
    da: spls-da analysis output
    treatments: considered treatments in the analysis
    datasetnames: the names of the input dataset
    annotations: annotations (if provided) of the elements in the analysis 
    danetwork: network obtained from splsda
}
\references{
mixOmics
}
\author{
Luis Valledor, Laura Lamelas
}
\seealso{
\code{\link{preprocess_omic_list}}
}
\examples{
testda<-da_analysis(listadatos = poldata,annotation = NULL,initialrow = 1,initialcolumn = 2,treatment1col = 1,treatment2col = NULL,treatment = 1,omiclevel = "proteome")}
