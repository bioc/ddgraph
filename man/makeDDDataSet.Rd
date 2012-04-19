\name{makeDDDataSet}
\alias{makeDDDataSet}
\title{Construct an DDDataSet object...}
\usage{makeDDDataSet(signal, name, classLabels, classLabelsCol, removeZeroVar=FALSE)
}
\description{Construct an DDDataSet object}
\value{a new DDDataSet object}
\arguments{\item{signal}{the matrix or data frame where rows are observations and columns variables}
\item{name}{the name of the dataset (to be used in plotting, etc)}
\item{classLabels}{the vector of class labels or target responses (aka target variable)}
\item{classLabelsCol}{the column which should be interpreted as class labels (either name or index)}
\item{removeZeroVar}{if to remove zero variance columns without producing an error (default: TRUE)}
}
\examples{# columns are features, rows observations
data <- matrix(rbinom(50, 1, 0.5), ncol=5)
# target class labels
labels <- c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1)
makeDDDataSet(data, name="example data", classLabels=labels)}

