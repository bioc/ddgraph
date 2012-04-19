\name{independent.contributions.formula.mul}
\alias{independent.contributions.formula.mul}
\title{Generate class labels by independent contributions of two variables}
\usage{independent.contributions.formula.mul(data, target.vars, prob1, prob0, logical.formula)
}
\description{Version of \code{independent.contributions.formula} that works with any number of variables. See the help page
for \code{independent.contributions.formula} for description of functionality.}
\value{a vector of binary class labels}
\arguments{\item{data}{a matrix or data.frame containing binary observations (columns are variables)}
\item{target.vars}{indexes of target variables}
\item{prob1}{vector of P(class labels = 1|varX=1) for different X}
\item{prob0}{vector of P(class labels = 1|varX=0) for different X}
\item{logical.formula}{a character string for the formula}
}
\examples{# noisy OR function with three variables and with noise level of 0.1 for a, b, and 0.2 for c
data <- cbind("a"=c(0,0,0,0,1,1,1,1), "b"=c(0,0,1,1,0,0,1,1), "c"=c(0,1,0,1,0,1,0,1))
independent.contributions.formula.mul(data, c("a", "b", "c"), c(0.9, 0.9, 0.8), c(0.1, 0.1, 0.2), "a | b | c")}

