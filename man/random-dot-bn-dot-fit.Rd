\name{random.bn.fit}
\alias{random.bn.fit}
\title{Generate a random \code{bn.fit} network}
\usage{random.bn.fit(nodes, num.neigh=2, prob.distr=prob.distr.norm, bn.graph)
}
\description{Generate a random Bayesian network using package \code{bnlearn}. The nodes specify the partial ordering
of the graph, and the conditional probabilities are sampled from given distribution. The network is
generated to have on average given number of neighbours (i.e. both in-going and out-going edges)}
\value{a list of two elements: \code{bn} - a \code{bn} object which contains the structure and 
\code{bn.fit} - a \code{bn.fit} object with filled in conditional probabilities}
\arguments{\item{nodes}{a vector of desired node names (basis for partial ordering)}
\item{num.neigh}{expected number of neighbours per node in the random graph}
\item{prob.distr}{the probability distribution function to use}
\item{bn.graph}{the \code{bn} object with an already laid out graph, if not supplied will be generated}
}
\examples{# a random network with 3 nodes "A", "B", "C" with average of 1 neighbour
random.bn.fit(c("A", "B", "C"), num.neigh=1)}

