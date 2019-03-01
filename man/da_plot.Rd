\name{da_plot}
\alias{da_plot}
\title{The function to plot the splsda analysis}
\description{
This function allows to fully depict the results from the spls-da analysis}
\usage{
da_plot()
}
\arguments{
  \item{da_list}{A "daanalysis" class object, see da_analysis function }
  \item{treatment}{numeric, c("1","2","3") dataset will be split by treatment 1; 2, dataset will be split by treatment 2; 3, dataset will be split by treatments1&2}
  \item{plottype}{c("Distance","Var","Composite","Composite_Annotated","Network")}
  \item{variablespace}{c("x","y","xy")}
  \item{useannot}{logical, if TRUE annotations will have}
  \item{fontsize}{fontsize in the plot and plot legend }
  \item{cutoff}{Network threshold value}
}
\details{
}
\author{
Luis Valledor}

\examples{
da_plot(da_list = daobject,treatment =1 ,plottype ="Distance" ,variablespace ="xy" ,useannot = FALSE ,cutoff = 0.85 )}
