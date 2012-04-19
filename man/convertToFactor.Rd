\name{convertToFactor}
\alias{convertToFactor}
\title{Convert data to factor representation}
\usage{convertToFactor(x)
}
\description{Convert a matrix, dataframe or vector into a factor representation. Each column is going
to be separately converted into a factor.}
\arguments{\item{x}{the input vector, data.frame or matrix}
}
\examples{# works on vectors, matrices and data frames
convertToFactor(0)
convertToFactor(c(1, 0, 0, 1, 0))
convertToFactor(matrix(c(1,0), nrow=2, ncol=2))
convertToFactor(data.frame("a"=c(1,0), "b"=c(0,1)))}

