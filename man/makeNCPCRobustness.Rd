\name{makeNCPCRobustness}
\alias{makeNCPCRobustness}
\title{Make a new NCPCRobustness object...}
\usage{makeNCPCRobustness(dataset, raw, params)
}
\description{Make a new NCPCRobustness object}
\details{Make a new NCPCRobustness object just with the raw resampling data and parameters used to generate them.
Should never directly use this function, but only via \code{DDDataSet::NCPCRobustness()}.}
\value{a new NCPCRobustness object}
\arguments{\item{dataset}{the DDDataSet object}
\item{raw}{the list of raw resampling classification of variables (direct, joint, etc..)}
\item{params}{the parameters used to generate the data (only the non-default one are listed)}
}

