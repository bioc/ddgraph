\name{toyExample}
\alias{toyExample}
\docType{data}

\title{
  A binary fictional toy example \code{DDDataSet} object.
}
\description{
  \code{toyExample} is an example dataset representing a set of 200 fictional cis-regulatory modules (CRMs). 
  The dataset contains binding patterns for two transcription factors A and B. It is used only in the package
  vignette.   
}

\usage{data(toyExample)}

\details{

In this fictional dataset we represent binding patterns of two transcription factors A and B on a set of CRMs. 
The target variable (T) is another binary vector that represents if a CRM is tissue specific or not (as obtained
by e.g. transgenic reporter assays). 

For more information and detailed examples see the package vignette. 
}

\examples{
data(toyExample)
calcDependence(toyExample)
}
\keyword{datasets}
