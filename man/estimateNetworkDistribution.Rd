\name{estimateNetworkDistribution}
\alias{estimateNetworkDistribution}
\title{Estimate network distribution parameters}
\usage{estimateNetworkDistribution(obj, use.class=FALSE)
}
\description{Estimate the in.degree distribution and conditional probability distribution from data}
\details{The algorithm uses hill-climbing with BIC to construct the network and estimate the
parameters. Then, provided that for each in-degree there is at least two nodes, it
estimates the beta distribution parameters.}
\value{a list of two elements: in.degree.distr - distribution of in-degrees, and 
beta.est - estimate beta distribution values}
\arguments{\item{obj}{and object of class DDDataSet}
\item{use.class}{if to include the class variable into the estimate}
}
\examples{data(mesoBin)
estimateNetworkDistribution(mesoBin$Meso)}

