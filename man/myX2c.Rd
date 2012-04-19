\name{myX2c}
\alias{myX2c}
\title{The Monte-Carlo chi-square test...}
\usage{myX2c(x, y, C, B=5000)
}
\description{The Monte-Carlo chi-square test}
\details{This is the reimplementation of Monte Carlo chi-square test to be sure
it works correctly. The Monte Carlo loop is implemented using \code{Rcpp} and uses
the R function \code{r2dtable()} to generate random contingency tables with
fixed marginals.}
\value{the P-value of the test}
\arguments{\item{x}{the first variable (vector of values)}
\item{y}{the second variable (vector of values)}
\item{C}{the variables to condition on - either a vector, or a list of vectors}
\item{B}{the number of Monte Carlo runs (defaults to 5000 if given NULL)}
}

