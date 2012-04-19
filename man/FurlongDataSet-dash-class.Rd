\name{FurlongDataSet-class}
\alias{FurlongDataSet-class}
\alias{FurlongDataSet}
\title{Data class for the Furlong dataset...}
\description{Data class for the Furlong dataset}
\details{A class to hold data from (Zizen 2009) paper (Supp Table 8). This class contains methods to convert it to
both binary and continuous \link{DDDataSet} objects.}
\references{Robert P. Zinzen et al., "Combinatorial binding predicts spatio-temporal cis-regulatory activity,"
Nature 462, no. 7269 (November 5, 2009): 65-70.}
\section{Slots}{\describe{\item{\code{signalMatrix}:}{(\code{\link{matrix}}) the signal matrix}\item{\code{targetClasses}:}{(\code{\link{factor}}) the target class names}}}
\section{Methods}{\describe{\item{\code{\link[=names,FurlongDataSet-method]{names}}}{\code{signature(x = "FurlongDataSet")}: Get the names of variables (column names of signal matrix)}\item{\code{\link[=signalMatrix,FurlongDataSet-method]{signalMatrix}}}{\code{signature(object = "FurlongDataSet")}: Retrieve the matrix with raw signal values}\item{\code{\link[=classLabels,FurlongDataSet-method]{classLabels}}}{\code{signature(object = "FurlongDataSet")}: Retrieve the vector of class labels (as factors)}\item{\code{\link[=toDDDataSet,FurlongDataSet-method]{toDDDataSet}}}{\code{signature(obj = "FurlongDataSet")}: Make the DDDataSet objects by selecting different tissues}}}

