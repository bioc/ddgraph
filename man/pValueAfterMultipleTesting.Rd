\name{pValueAfterMultipleTesting}
\alias{pValueAfterMultipleTesting}
\title{Multiple testing correction procedure for ncpc()}
\usage{pValueAfterMultipleTesting(dsep, x, adjC.pvals.at.n, p.value.adjust.method)
}
\description{This function is only for DDGraph with multiple testing correction enabled. The overall procedure is similar to
that described in (Li&Wang 2009). This is a helper function for \code{DDDataSet:ncpc()}. The single P-value
of D-separation is substituted in the list of P-values, P-values adjusted and the resulting P-value after correction
in the context of other P-values reported.}
\value{the p value after multiple test correction (if any)}
\references{J. Li and Z. J Wang, "Controlling the false discovery rate of the association/causality structure 
learned with the PC algorithm" The Journal of Machine Learning Research 10 (2009): 475-514.}
\arguments{\item{dsep}{the conditional independence test result (of type \code{CITestResult})}
\item{x}{the index of the variables}
\item{adjC.pvals.at.n}{the p values associated with the variables at size n of conditioning set (list [[n]] -> [pvals])}
\item{p.value.adjust.method}{the p value adjustment method (same as in p.adjust())}
}

