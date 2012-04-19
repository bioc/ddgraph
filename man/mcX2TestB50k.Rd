\name{mcX2TestB50k}
\alias{mcX2TestB50k}
\title{Wrapper around the bnlearn mc-x2 test (B=50k)}
\usage{mcX2TestB50k(x, y, S, suffStat)
}
\description{Version of \code{mcX2Test()} with 50000 Monte Carlo replicates.}
\value{p value of the test}
\arguments{\item{x}{the index of the first variable}
\item{y}{the index of the second variable}
\item{S}{the conditioning set}
\item{suffStat}{the sufficient statistics to do the test, in this case a list of one element: 
dm where the values matrix is stored}
}
\examples{suffStat <- list(dm = cbind("a"=c(0,1,0,0,1,0), "b"=c(1,0,0,0,1,0), "c"=c(0,0,0,1,1,1)))
# test if a is independent of b
mcX2TestB50k(1, 2, NULL, suffStat)
# test if a is independent of b conditioned on c
mcX2TestB50k(1, 2, 3, suffStat)}

