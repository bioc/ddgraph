\name{biased.bn.fit}
\alias{biased.bn.fit}
\title{Random network with a biased degree distribution}
\usage{biased.bn.fit(nodes, beta.est, in.degree.distr, bn.graph)
}
\description{A version of random.bn.fit which generates a graph based on degree distribution and beta distribution for probabilities}
\value{a list of two elements: \code{bn} - a \code{bn} object which contains the structure and 
\code{bn.fit} - a \code{bn.fit} object with filled in conditional probabilities}
\arguments{\item{nodes}{character vector of node names}
\item{beta.est}{the beta distribution parameters for different degrees of a node. Should be a list where [[2]] corresponds 
to 2-dimenstional contingency table (i.e. one parent, one output). It contains a data.frame with columns shape1, shape2
for the beta distribution, and rows are degrees of freedom (in this case 2, when P(Out=0|Parent=0) and P(Out=0|Parent=1))}
\item{in.degree.distr}{a vector with degree distribution for all the nodes in the network (names are ignored, and degree is 
randomly sampled from this vector)}
\item{bn.graph}{if the graph structure is already available, then the graph structure in object of class "bn"}
}
\examples{# nodes, conditional probability distribution, an indegree distribution
nodes = letters[1:5]
beta.est = list(data.frame(shape1=2,shape2=3), data.frame(shape1=c(2,4), shape2=c(5,2)), data.frame(shape1=c(1,2,3,4), shape2=c(3,2,1,2)))
in.degree.distr = c(0, 1, 1, 2, 2)
# make a random graph using these parameters
biased.bn.fit(nodes, beta.est, in.degree.distr)}

