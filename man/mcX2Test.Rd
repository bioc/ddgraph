\name{mcX2Test}
\alias{mcX2Test}
\title{Wrapper around the bnlearn mc-x2 test}
\usage{mcX2Test(x, y, S, suffStat)
}
\description{Implements the mc-x2 test in format needed for pcalg.}
\value{p value of the test}
\arguments{\item{x}{the index of the first variable}
\item{y}{the index of the second variable}
\item{S}{the conditioning set}
\item{suffStat}{the sufficient statistics to do the test, in this case a list of one element: 
dm where the values matrix is stored}
}
\examples{suffStat <- list(dm = cbind("a"=c(0,1,0,0,1,0), "b"=c(1,0,0,0,1,0), "c"=c(0,0,0,1,1,1)))
# test if a is independent of b
mcX2Test(1, 2, NULL, suffStat)
# test if a is independent of b conditioned on c
mcX2Test(1, 2, 3, suffStat)}

