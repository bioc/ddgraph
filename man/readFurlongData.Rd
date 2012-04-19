\name{readFurlongData}
\alias{readFurlongData}
\title{Read the Furlong Dataset}
\usage{readFurlongData(infile)
}
\description{Read the Furlong data into a FurlongDataSet object.}
\details{Read the Furlong Dataset form the Supplementary Table 8 file provided with the package. An alternative
filename can be specified as well.}
\value{an object of type FurlongDataSet witht the loaded data}
\arguments{\item{infile}{the filename to load from, default to supplementary_table_8_training_set.txt in extdata/ of package}
}
\examples{# read the furlong dataset that is provided with the package
readFurlongData()}

