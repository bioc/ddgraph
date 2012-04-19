\name{extractCITestResultProperty}
\alias{extractCITestResultProperty}
\title{Extract CITestResult properties}
\usage{extractCITestResultProperty(ciTestList, prop.name)
}
\description{This is a helper function for \code{DDDataSet::ncpc()}. From a list of ciTestResult object extract 
a list containing only one property.}
\value{a vector with the extracted property}
\arguments{\item{ciTestList}{a two-level list of ciTestResult objects}
\item{prop.name}{the name of the property to extract (one of the slot names)}
}

