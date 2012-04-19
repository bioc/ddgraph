\name{CITestResult-class}
\alias{CITestResult-class}
\alias{CITestResult}
\title{Data class to store the results of a conditional independence test...}
\description{Data class to store the results of a conditional independence test}
\details{This class stored the results from \code{DDDataSet::ciTest()}. It stores the
indexes and names of two variables involved in the test, the conditioning set
as well as the P-value and type of test.}
\section{Slots}{\describe{\item{\code{targetInx}:}{(\code{\link{numeric}}) the index of the first variable}\item{\code{targetName}:}{(\code{\link{character}}) the name of the first variable}\item{\code{sourceInx}:}{(\code{\link{numeric}}) the index of the second variable}\item{\code{sourceName}:}{(\code{\link{character}}) the name of the second variable}\item{\code{condSetInx}:}{(\code{\link{numeric}}) the indexes of variables we condition on}\item{\code{condSetName}:}{(\code{\link{character}}) the names of variables we condition on}\item{\code{pValue}:}{(\code{\link{numeric}}) the associated p value}\item{\code{testType}:}{(\code{\link{character}}) the type of the conditional independence test performed}\item{\code{reliable}:}{(\code{\link{logical}}) if this appears to be a reliable test of conditional independence}}}
\section{Methods}{\describe{\item{\code{\link[=$,CITestResult-method]{$}}}{\code{signature(x = "CITestResult")}: Access slots using the dollar notation}\item{\code{\link[=[[,CITestResult,ANY,ANY-method]{[[}}}{\code{signature(x = "CITestResult", i = "ANY", j = "ANY")}: Access slots using the double square bracket notation}\item{\code{\link[=names,CITestResult-method]{names}}}{\code{signature(x = "CITestResult")}: Names of slots that can be accessed with $ notation}\item{\code{\link[=show,CITestResult-method]{show}}}{\code{signature(object = "CITestResult")}: show method for CITestResult}}}

