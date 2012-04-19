\name{formulaFalseNeg}
\alias{formulaFalseNeg}
\title{Generate class labels by a noisy formula with high false negative rate}
\usage{formulaFalseNeg(data, var1, var2, false.neg, logical.formula)
}
\description{Generate class labels by using the readout mechanism. Logical formula is applied to two variables
which are read out from the real data using the var1 and var2 probabilities. This only works
with binary variables.}
\value{a binary vector containing the class labels}
\arguments{\item{data}{a matrix or data.frame containing binary observations (columns are variables)}
\item{var1}{index or name of the first variable}
\item{var2}{index or name of the second variable}
\item{false.neg}{a false negative probability}
\item{logical.formula}{logical formula to apply}
}
\examples{# noisy OR function with 0.1 probability of error for reading "a" and "b" (error in both 1 and 0)
data <- cbind("a"=c(0,0,1,1), "b"=c(0,1,0,1))
formulaFalseNeg(data, "a", "b", 0.8, "a | b")}

