% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SPADEVizR_unimodality.R
\name{SPADEVizR_unimodality}
\alias{SPADEVizR_unimodality}
\title{#########################################################################################
Test unimodality of metaclusters using SPADEVizR function}
\source{
https://academic.oup.com/bioinformatics/article/33/5/779/2662307

https://github.com/tchitchek-lab/SPADEVizR
}
\usage{
SPADEVizR_unimodality(df = NULL, cluster = NULL, clusters = NULL,
  clustering.markers = NULL, uniform.test = "unimodality",
  th.pvalue = 0.05, th.IQR = 2, density.PDFfile = NULL,
  density.PDFfile.dim = c(8, 5), heatmap.PDFfile = NULL,
  tile.color = "black", verbose = TRUE)
}
\arguments{
\item{df}{a matrix of numerical values}

\item{cluster}{a null value}

\item{clusters}{a null value}

\item{clustering.markers}{an array of character string}

\item{uniform.test}{a character string}

\item{th.pvalue}{a numerical value}

\item{density.PDFfile}{a character string}

\item{density.PDFfile.dim}{an array of 2 numerical values}

\item{heatmap.PDFfile}{a character string}

\item{verbose}{a logical}
}
\value{
uniform.res
}
\description{
Patrice Vallin, easyFlow, Sept 2019.
}
\details{
Adapted from :
Guillaume Gautreau et al., Bioinformatics, March 2017
https://academic.oup.com/bioinformatics/article/33/5/779/2662307
https://github.com/tchitchek-lab/SPADEVizR
########################################################################################
}
