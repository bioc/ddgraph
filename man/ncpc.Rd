\name{ncpc}
\alias{ncpc}
\title{Make a Direct Dependence Graph using the NCPC algorithm...}
\usage{ncpc(obj, alpha=0.05, p.value.adjust.method="none", test.type=c("mc-x2-c",
    "cor"), max.set.size=NULL, mc.replicates=5000, report.file=NULL,
    verbose=FALSE, star=FALSE, min.table.size=10)
}
\description{Make a Direct Dependence Graph using the NCPC algorithm}
\details{Make a Direct Dependence Graph using a P-value and conditional independence tests. There are two version of the
algorithm: NCPC and NCPC*. NCPC finds the causal neighbourhood while the NCPC* infers the full Markov Blanket. 

The full algorithm is given in (Stojnic et al, 2012).}
\value{DDGraph object}
\references{R. Stojnic et al (2012): "A Graphical Modelling Approach to the Dissection of Highly Correlated Transcription Factor Binding Site Profiles", in press, PloS Computational Biology.}
\arguments{\item{obj}{DDDataSet object}
\item{alpha}{the alpha (P-value) cutoff for conditional independence tests}
\item{p.value.adjust.method}{the multiple testing correction adjustment method}
\item{test.type}{the type of conditional independence test (default: Monte Carlo x2 test "mc-x2-c" for binary data
and partial correlation "cor" for continuous data) . See the documentation 
for \code{\link{ciTest}} for other available conditional independence tests}
\item{max.set.size}{the maximal number of variables to condition on, if NULL
estimated from number of positives in class labels (default: NULL)}
\item{mc.replicates}{the number of Monte-carlo replicates, if applicable (default: 5000)}
\item{report.file}{name of the file where a detailed report is to be printed, 
reporting is suppressed if NULL (default: NULL)}
\item{verbose}{if to print out information about how the algorithm is progressing (default: TRUE)}
\item{star}{if to use the NCPC* algorithm (default: FALSE)}
\item{min.table.size}{the minimal number of samples in a contingency table per conditioning set 
(makes sense only for discrete data)}
}
\examples{### load binary data for Mesoderm
data(mesoBin)
# run the NCPC algorithm with alpha=0.05 (on discrete data)
ncpc(mesoBin$Meso, alpha=0.05, test.type="mc-x2-c")
# run the NCPC* algorithm with alpha=0.05 (on discrete data)
res <- ncpc(mesoBin$Meso, alpha=0.05, test.type="mc-x2-c", star=TRUE)

# analysis of results:
class(res)
# although of class DDGraph, behaves much like a list
names(res)
# parameters used in obtaining results
res$params
# labels for each of the variables
res$final.calls
# direct variables
res$direct

### load continous data
data(mesoCont)
# run the NCPC algorith with alpha=0.05 (on continuous data)
ncpc(mesoCont$Meso, alpha=0.05, test.type="cor", max.set.size=1)
# run the NCPC* algorith with alpha=0.05 (on continuous data)
ncpc(mesoCont$Meso, alpha=0.05, test.type="cor", max.set.size=1, star=TRUE)}

