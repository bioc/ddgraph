\name{is.binary}
\alias{is.binary}
\title{Check if data structure has binary data in it}
\usage{is.binary(x)
}
\description{Check if a vector, data frame or matrix contains only binary (0,1) values.}
\value{boolean TRUE or FALSE}
\arguments{\item{x}{the input vector, data.frame or matrix}
}
\examples{# works on vectors, matrices and data frames
is.binary(0)
is.binary(c(1, 0, 0, 1, 0))
is.binary(matrix(c(1,0), nrow=2, ncol=2))
is.binary(data.frame("a"=c(1,0), "b"=c(0,1)))

# returns FALSE if not binary
is.binary(c(1, 2, 3))}

