% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mapman_plot.R
\name{mapman_plot}
\alias{mapman_plot}
\title{Generate plots based on Mapman categories and groups (mapman_group ())}
\usage{
mapman_plot(MMOlist, plottype = "bar", ploterror = TRUE, logscale = FALSE,
hmparameters = c("manhattan","manhattan","ward.D",TRUE,"row"), normbyelements = FALSE,
fontsizes = c(14,10,16,12))
}
\arguments{
\item{MMOlist}{Mapman aggrupation table with sums of abundances according to mapman categories}

\item{plottype}{Character.Type of plot to be displayed, "bar", "heatmap" or "pie".}

\item{ploterror}{Boolean indicating whether errors within replicas should be shown}

\item{logscale}{Boolean indicating whether axis of bar plot should be log-scaled}

\item{hmparameters}{Character vector (length 5) determining main features of heatmap plot (if defined as plottype argument).
hmparameters <- c("distance measure used in row clustering","distance measure in column clustering","clustering method","displaynumbers","scale")
where distance measure used in row/col clustering: character, those methods supported by stats::dist ("correlation","euclidean","maximum","manhattan","camberra","binary","minkowski")
clustering method: character, those methods provided by stats::hclust ("ward.D","ward.D2","single","complete","average","mcquiry","median","centroid")
displaynumbers: boolean determining if numeric values should be indicated in heatmap cells.
scale: character, direction of abundance balancing ("row", "column","none"). If set to row or col the given numbers represent the percentage by bin or by sample respectively.}

\item{normbyelements}{Boolean indicating if MapMan classification abundances should be divided by the number of elements of each bin. Default FALSE}

\item{fontsizes}{Numeric vector, (length = 4) containing fontsizes of title, subtitle and axis labels}
}
\value{
Plot
}
\description{
A function which is used to plot abundances of mapman_group object. Barplot, pie charts and heatmaps available
}
\note{
Mapman annotation files could contain duplicates, importannotation function will ask you how to proceed.
}
\author{
Luis Valledor and Laura Lamelas
}
