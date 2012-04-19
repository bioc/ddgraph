\name{toDDDataSet,FurlongDataSet-method}
\alias{toDDDataSet,FurlongDataSet-method}
\alias{toDDDataSet}
\title{DDDataSet object from FurlongDataSet}
\usage{\S4method{toDDDataSet}{FurlongDataSet}(obj, tissues=c(), convertToBinary=TRUE, prettyNames=FALSE, ...)
}
\description{Make the DDDataSet objects by selecting different tissues}
\value{either single DDDataSet object, or a list of them (depending on number of selected tissues)}
\arguments{\item{obj}{the FurlongDataSet object}
\item{tissues}{tissue names for which DDDataSet objects should be generated (default to all available tissues)}
\item{convertToBinary}{if to convert the signal into binary values}
\item{prettyNames}{if to make the names pretty, e.g. twi_2.4 -> Twi 2-4h}
\item{...}{unused}
}
\examples{# load binarized data with prettified names
all.data <- toDDDataSet(readFurlongData(), prettyNames=TRUE)
# load continuous data with original names
all.data <- toDDDataSet(readFurlongData(), convertToBinary=FALSE)}

