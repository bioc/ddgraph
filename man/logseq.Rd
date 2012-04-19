\name{logseq}
\alias{logseq}
\title{Generate sequence in log scale}
\usage{logseq(a, b, n=8)
}
\description{Generate sequence but in log scale. This function takes takes the length of log-sequence and the minimal and
maximal point. It returns the interval between \code{a} and \code{b} divided in log scale.}
\value{a vector of numbers}
\arguments{\item{a}{the smaller value in the interval}
\item{b}{the bigger value in the interval}
\item{n}{the number of intervals to divide a,b into}
}
\examples{# produces vector c(0.01, 0.1, 1)
logseq(0.01, 1, 3)}

