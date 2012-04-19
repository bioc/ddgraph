\name{blockingNodes}
\alias{blockingNodes}
\title{Find all such nodes in neighbourhood of source node that are blocking at least one active path leading to another node...}
\usage{blockingNodes(allPaths, nodes)
}
\description{Find all such nodes in neighbourhood of source node that are blocking at least one active path leading to another node}
\value{a list with blocking nodes and minimal length to the target node: target node => blocked by => number of steps}
\arguments{\item{allPaths}{a list of active paths from a source node (as produced by activePaths())}
\item{nodes}{a vector of target nodes for which we are finding blocking nodes}
}

