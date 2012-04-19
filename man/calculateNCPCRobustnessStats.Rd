\name{calculateNCPCRobustnessStats}
\alias{calculateNCPCRobustnessStats}
\title{Calculate NCPCRobustness statistics...}
\usage{calculateNCPCRobustnessStats(obj)
}
\description{Calculate NCPCRobustness statistics}
\details{Calculate the statistics for the NCPCRobustness object - this is separate from object 
construction for convenience of testing, should always be called after object creation. 
Never use directly (except for testing), use instead via \code{DDDataSet::NCPCRobustness()}.}
\value{the modified NCPCRobustness object with the statistics calculated}
\arguments{\item{obj}{NCPCRobustness object}
}

