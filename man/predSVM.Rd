\name{predSVM}
\alias{predSVM}
\title{Calculate the decision value of an SVM model}
\usage{predSVM(f, feature)
}
\description{Calculate the decision value of an SVM model. Note this is different from the actual
prediction which is either 0 or 1, while decision values go from -1 to 1.
(taken from [Zizen 2009] supplementary code)}
\value{Decision value in the range -1 to 1.}
\arguments{\item{f}{The trained SVM model object.}
\item{feature}{The input value to which output is needed.}
}

