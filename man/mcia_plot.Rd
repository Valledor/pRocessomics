% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mcia_plot.R
\name{mcia_plot}
\alias{mcia_plot}
\title{Multiple coinertia plot}
\usage{
mcia_plot(mcia_list, treatment = 1, compX = 1,
compY = 2, compZ = NULL, plottype = "Composite", useannot = FALSE,
levelchoice = "all", fontsizes = c(14,10,16,12), fortopscoring = c(1,20,"abs"),
confidence = 0.9)
}
\arguments{
\item{mcia_list}{List containing the multiple coinertia analysis, generated with pRocessomics::mcia_analysis function}

\item{treatment}{numeric value to determine the colors and symbols to be displayed in the plot, used in Score plot and Topscring plot. Available optiones are:
\itemize{
\item 1 to show different colors according to treatment 1 (treatment1col)
\item 2 to show different colors according to treatment 2 (treatment2col)
\item 3 to show different colors according to the combination of both treatments
\item 4 to show different colors according to treatment 1 and different symbols for treatment 2
}}

\item{compX}{component to be plotted in X axis in Score and Var plots, default 1}

\item{compY}{component to be plotted in X axis in Score and Var plots, default 2}

\item{compZ}{component to be plotted in Z axis for 3D plots in SynVar and Var plots, default NULL}

\item{plottype}{available options:
\itemize{
\item Scree, explained variance plot
\item SynVar, sample plot each dot represent a sample
\item SynVar_Ellipse, score plot with confidence ellipses by treatment
\item Topscoring, topscoring variables in the selected dimension
\item Var, variable plot, each dot represent a variable
\item Biplot, biplot showing samples and variables
\item CoIn, coinertia plot showing samples and omiclevels
\item Pseudoeigen, plot showing correlations among datasets
\item Composite: combination of the four most representative plots of the analysis
\item All all aforementioned plots will be exported as .pdf files in your working directory
}}

\item{useannot}{logical, indicating whether an annotation file should be used. Please note that the annotation file can be only generated by pRocessomics::importannotation() function}

\item{levelchoice}{character vector; quoted name(s) of omiclevels (one or more) to plot, can be set to "all" for include every omiclevel}

\item{fontsizes}{vector containing the font sizes to use in the plots}

\item{fortopscoring}{vector (length = 3) indicating the number of the component to be analyzed, the top n loadings to plot and the mode to compute "abs" for absolute value and "posneg" for positive and negative loading analyzed separately}

\item{confidence}{numeric value between 0 and 1 for confidence ellipses in SynVar_Ellipse plot}
}
\value{
mcia plot
}
\description{
A function to plot different aspects of multiple coinertia analysis
}
\author{
Luis Valledor and Laura Lamelas
}
