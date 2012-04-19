\name{DDGraphEdge-class}
\alias{DDGraphEdge-class}
\alias{DDGraphEdge}
\title{An edge in an DDGraph...}
\description{An edge in an DDGraph}
\details{This class represents an edge in an Direct Dependence Graph. It is normally found in the
\code{DDGraph::edges} list. It records the source and target nodes for the edge, the
edge type, as well as the conditional independence tests it represents.}
\section{Slots}{\describe{\item{\code{fromInx}:}{(\code{\link{numeric}}) the index of the first variable from which the edge goes}\item{\code{fromName}:}{(\code{\link{character}}) the name of the first variable from which the edge goes}\item{\code{toInx}:}{(\code{\link{numeric}}) the index of the second variable to which the edge goes}\item{\code{toName}:}{(\code{\link{character}}) the name of the second variable to which the edge goes}\item{\code{ciTests}:}{(\code{\link{list}}) a list of associated \code{CITestResult} objects}\item{\code{type}:}{(\code{\link{character}}) type of edge: "directed", "undirected", "bidirectional", "dashed"}}}
\section{Methods}{\describe{\item{\code{\link[=show,DDGraphEdge-method]{show}}}{\code{signature(object = "DDGraphEdge")}: show method for DDGraphEdge}}}

