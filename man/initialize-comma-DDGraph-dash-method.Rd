\name{initialize,DDGraph-method}
\alias{initialize,DDGraph-method}
\title{Construct new DDGraph object...}
\usage{\S4method{initialize}{DDGraph}(.Object, ..., direct=vector(mode = "numeric"), indirect=vector(mode =
    "numeric"), joint=vector(mode = "numeric"), conditional=vector(mode
    = "numeric"), conditionalJoint=vector(mode = "numeric"),
    edges=list(), dataset=new("DDDataSet"), params=list(), stats=list())
}
\description{Construct new DDGraph object}
\details{Properly initialize the object}
\arguments{\item{.Object}{DDGraph object}
\item{direct}{direct variable indexes}
\item{indirect}{indirect variable indexes}
\item{joint}{joint variable indexes}
\item{conditional}{conditional variable indexes}
\item{conditionalJoint}{conditionally joint variable indexes}
\item{edges}{edges list}
\item{dataset}{DDDataSet object}
\item{params}{parameters used to make this object}
\item{stats}{the statistics used to make this object}
\item{...}{unused}
}

