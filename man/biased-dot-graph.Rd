\name{biased.graph}
\alias{biased.graph}
\title{Generate random network with degree distribution}
\usage{biased.graph(nodes, in.degree.distr)
}
\description{Generate a random directed graph with the given node ordering and degree distribution}
\value{an object of class bn with the random graph}
\arguments{\item{nodes}{character vector of node names which species the node ordering}
\item{in.degree.distr}{the node in-degree distribution}
}
\examples{# a random network of 5 nodes with choosen in-degree distribution
biased.graph(letters[1:5], c(0, 1, 1, 2, 2))}

