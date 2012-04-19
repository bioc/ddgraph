\name{prob.distr.norm}
\alias{prob.distr.norm}
\title{Normal distribution function for \code{random.bn.fit}}
\usage{prob.distr.norm(n, sd=1/3)
}
\description{Generate 2^n numbers from distribution with most of the pdf mass in extreme probabilities (mirrored normal). 
We use standard deviation of 1/3 and modulo-1 of normal distribution.}
\value{vector of 2^n random numbers}
\arguments{\item{n}{number of variables}
\item{sd}{the standard deviation of distribution}
}
\examples{# return 8 random numbers since n=3
prob.distr.norm(3)}

