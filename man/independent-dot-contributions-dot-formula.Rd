\name{independent.contributions.formula}
\alias{independent.contributions.formula}
\title{Generate class labels by independent contributions of two variables}
\usage{independent.contributions.formula(data, var1, var2, var1.prob1, var1.prob0, var2.prob1, var2.prob0,
    logical.formula, false.neg=0, false.pos=0)
}
\description{Generate class labels by using the readout mechanism. Logical formula is applied to two variables
which are read out from the real data using the var1 and var2 probabilities. This only works
with binary variables.}
\value{a binary vector containing the class labels}
\arguments{\item{data}{a matrix or data.frame containing binary observations (columns are variables)}
\item{var1}{index or name of the first variable}
\item{var2}{index or name of the second variable}
\item{var1.prob1}{the conditional probability P(class labels = 1|var1=1)}
\item{var1.prob0}{the conditional probability P(class labels = 1|var1=0)}
\item{var2.prob1}{the conditional probability P(class labels = 1|var2=1)}
\item{var2.prob0}{the conditional probability P(class labels = 1|var2=0)}
\item{logical.formula}{logical formula to apply}
\item{false.neg}{a false negative probability}
\item{false.pos}{a false positive probability}
}
\examples{# noisy OR function with 0.1 probability of error for reading "a" and "b" (error in both 1 and 0)
data <- cbind("a"=c(0,0,1,1), "b"=c(0,1,0,1))
independent.contributions.formula(data, "a", "b", 0.9, 0.1, 0.9, 0.1, "a | b")}

