\name{combinationsTest}
\alias{combinationsTest}
\title{Significant combinations of variables}
\usage{combinationsTest(obj, selected.vars, cutoff=0.05, p.adjust.method="none", verbose=TRUE)
}
\description{Calculate which combinations of values of variables are significantly different in the 
two classes (only for binary data). This function takes an DDDataSet and a number of variables
and finds those combinations of values of those variables that have significantly different
frequencies in the two class labels.}
\value{\code{data.frame} with ordered combinatorial patterns of selected variables}
\arguments{\item{obj}{DDDataSet object}
\item{selected.vars}{indexes or names of variables selected for the test}
\item{cutoff}{the p value cutoff for reporting (default: 0.05)}
\item{p.adjust.method}{the multiple adjustment method (default: none)}
\item{verbose}{if to print progress output and additional information}
}
\examples{data(mesoBin)
# find significant differences at 0.2 FDR
combinationsTest(mesoBin$Meso, c("Twi 2-4h", "Tin 6-8h", "Mef2 6-8h"), 0.2, "fdr")}

