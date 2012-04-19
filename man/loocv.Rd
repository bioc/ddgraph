\name{loocv}
\alias{loocv}
\title{Leave-one-out cross validation}
\usage{loocv(data, train.fun, eval.fun, verbose=FALSE)
}
\description{Leave-one-out cross validation systematically leaves out one row from the data, retrains the
classifier and then uses the retrained classifier to make a prediction for the left-out row.}
\value{A vector of length \code{nrow(data)} containing predictions from \code{eval.fun} when each 
row is left out once}
\arguments{\item{data}{The \code{data.frame} with data. Columns are variables, rows are observations.}
\item{train.fun}{The training function that takes the data without one of the rows left out.}
\item{eval.fun}{The prediction function that takes the trained model and the left out data point.}
\item{verbose}{If to print progress indication}
}

