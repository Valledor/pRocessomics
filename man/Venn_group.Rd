\name{Venn_group}
\alias{Venn_group}
\title{Venn diagram Table in Omic Datasets}

\description{
This function processes a POL object with the different omic levels to be analized in order to provide the Venn diagram's table.
}

\usage{
Venn_group(datalist,
                    initialrow=1,
                    initialcolumn=3,
                    treatment1col=1,
                    treatment2col=2,
                    treatment=1,
                    omiclevel=NULL,
                    threshold=NULL)
                    }
\arguments{
\item{datalist}{A "POL" object.}
\item{initialrow }{First row of numeric data within datasets.}
\item{initialcolumn}{First column of numeric data within datasets.}
\item{treatment1col}{Column number in which treatment 1 is indicated.}
\item{treatment2col}{Column number in which treatment 2 is indicated (if applicable). If only one treatment is set indicate the same column number than treatment1col}
\item{treatment}{1 dataset will be split by treatment 1; 2, dataset will be split by treatment 2; 3, dataset will be split by treatments1&2}
\item{omiclevel}{Concrete omic layer to make the anaylisis. If set to NULL or "all" all data in datalist will be analysed together}
\item{threshold}{Minimum number of not-null values required to take into account one variable  }
}

\value{
  Returns a "venannalysis" object, a list consisting of two dataframes, "treatments" and "data", which contain respectively the treatments present in the initial datalist and the final Venn diagram table with the Venn region, the frecuency, and the variables present in that region.
}


\details{
The aim of \code{Venn_group} is to provide a \code{vennanalysis} object containing Venn regions and the list of variables in each region, required to make a Venn diagram plot. This function requires a \code{POL} object, previously defined by using \code{preprocess_omic_list}.
}


\author{
  Luis Valledor, Laura Lamelas
}

\references{
  \code{\link{citasd}}.
}


\seealso{
\code{\link{Venn_plot}},
\code{\link{preprocess_omic_list}},
and http://www.miweb.es for more details.
}

\examples{
data(xxx)
my.Venn.table <- Venn_group (my_preprocessed_list, initialrow = 1, initialcolumn = 3, treatment1col = 1, treatment2col = 2, treatment = 3, omiclevel="Proteins", threshold=1)
}
