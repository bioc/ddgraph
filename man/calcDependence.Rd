\name{calcDependence}
\alias{calcDependence}
\title{Dependence with target variable}
\usage{calcDependence(dd, method="ncpc", ...)
}
\description{Calculate dependence with a target variable}
\details{This function is a front-end convenience function to access predictions of direct
dependence with a target variable by various Graphical Modelling algorithm.  

Consider a set of variable X_1, ..., X_m and a target variable T. We say that
that X_i is directly dependent with T if there is no other set of variable X_j, X_k, ...
such that it renders X_i conditionally independent of T. In other words, 
X_i is the most immediate casual cause/consequence of T in the set of chosen variables. 

Note that the above statement is different from that of classical feature selection for
classification. A set of features obtained with feature selection have the property that
a good classifier can be made based on them alone, while the above statement establishes 
statistical properties of variables. The set of variables with direct dependence
might not be optimal for classification, since classification performance can be
strongly influenced by false negatives (Friedman et al, 1997).}
\references{Nir Friedman, Dan Geiger, and Moises Goldszmidt, "Bayesian Network Classifiers", 
Machine Learning 29 (November 1997): 131-163.}
\value{A list with elements:
\itemize{
\item obj - the resulting object, either of class \code{DDGraph} for ncpc and ncpc* algorithms, or
of class \code{bn} for \code{bnlearn} algorithms, or 
of class \code{pcAlgo} for PC algorithm.
\item nbr - the variables with direct dependence (i.e. target node neighbourhood in the causal graph). 
For both ncpc and ncpc* includes variables with direct and joint dependence. 
\item mb - the variables in Markov Blanket of target variable. Not applicable for ncpc algoritm. For ncpc* 
algorithm includes variables with direct, joint and conditional dependence.
\item labels - for ncpc and ncpc* contains the set of labels that are output of the algorithm. 
}}
\arguments{\item{dd}{An object of type DDDataSet}
\item{method}{Algorithm to use. Valid values are:
\itemize{
\item ncpc - Neighbourhood Consistent PC algorithm
\item ncpc* - Neighbourhood Consistent PC algorithm star version
\item hc - Hill-climbing with custom penalty functions
\item hc-bic - Hill-climbing with BIC penalization (package bnlearn)
\item hc-bde - Hill-climbing with BDe penalization (package bnlearn)
\item iamb - IAMB algorithm (package bnlearn)
\item fast.iamb - FastIAMB algorithm (package bnlearn)
\item inter.iamb - InterIAMB algorithm (package bnlearn)
\item pc - PC algorithm (package pcalg)
\item mmpc - MMPC algorithm (package bnlearn)
\item mmhc - MMHC with custom penalty functions
\item mmhc-bic - MMHC with BIC penalization (package bnlearn)
\item mmhc-bde - MMHC with BDe penalization (package bnlearn)
}}
\item{...}{Extra parameters passed to backend functions \code{ncpc()}, \code{plotBNLearn()} 
and \code{plotPCAlgo()} depending on the picked algorithm (parameter \code{method}).

Extra parameters for ncpc and ncpc*:
\itemize{
\item alpha - the alpha (P-value) cutoff for conditional independence tests (default: 0.05)
\item p.value.adjust.method - the multiple testing correction adjustment method (default: "none")
\item test.type - the type of conditional independence test (default: "mc-x2-c"). See the documentation 
for \code{\link{ciTest}} for available conditional independence tests
\item max.set.size - the maximal number of variables to condition on, if NULL
estimated from number of positives in class labels. Needs to be specified for 
continuous data. (default: NULL)
\item mc.replicates - the number of Monte-Carlo replicates for the conditional independence
test, if applicable (default: 5000)
\item report.file - name of the file where a detailed report is to be printed, 
reporting is suppressed if NULL (default: NULL)
\item verbose - if to print out information about how the algorithm is progressing (default: TRUE)
\item min.table.size - the minimal number of samples in a contingency table per conditioning set 
(applicable only for discrete data) (default: 10)
}  
Extra parameters for hc, mmhc:
\itemize{
\item score - score function to use, accepts all from bnlearn package. For discrete 
data: "loglik", "aic", "bic", "bde", "k2". For continuous: "loglik-g", "aic-g", "bic-g",
"bge". For more details see help page for \code{package-bnlearn}. 
\item make.plot - if to make a plot or just return the network (default: FALSE)
\item blacklist - a data frame with two columns (optionally labeled "from" and "to"), 
containing a set of arcs not to be included in the graph. (default: NULL)
\item restart - the number of random restarts for score-based algorithms (default: 0)
\item scale - the colour scaling (default: 1.5)
\item class.label - the label to use for the target variable (default: "target")
\item use.colors - if to colour code the enrichment/depletion in a plot (default: TRUE)
}
Extra parameters for hc-bic, hc-bde, mmhc-bic, mmhc-bde:
\itemize{
\item make.plot - if to make a plot or just return the network (default: FALSE)
\item blacklist - a data frame with two columns (optionally labeled "from" and "to"), 
containing a set of arcs not to be included in the graph. (default: NULL)
\item restart - the number of random restarts for score-based algorithms (default: 0)
\item scale - the colour scaling (default: 1.5)
\item class.label - the label to use for the target variable (default: "target")
\item use.colors - if to colour code the enrichment/depletion in a plot (default: TRUE)
}
Extra parameters for iamb, fast.iamb, inter.iamb, mmpc:
\itemize{
\item make.plot - if to make a plot or just return the network (default: FALSE)
\item alpha - the alpha value of conditional independence tests (default: 0.05)
\item test - the type of conditional independence test (default: "mc-mi"). For 
conditional independence tests available consult the bnlearn package help 
page (?bnlearn). 
\item B - the number of Monte-Carlo runs for conditional independence tests,
if applicable (default: 5000)
\item blacklist - a data frame with two columns (optionally labeled "from" and "to"), 
containing a set of arcs not to be included in the graph. (default: NULL)
\item scale - the colour scaling (default: 1.5)
\item class.label - the label to use for the target variable (default: "target")
\item use.colors - if to colour code the enrichment/depletion in a plot (default: TRUE)
}
Extra parameters for pc:
\itemize{
\item alpha - the alpha value cut-off for the conditional independence tests (default: 0.05)
\item verbose - if to show progress (default: FALSE)
\item directed - if TRUE applies PC algorithm, if FALSE applies PC-skeleton (default: TRUE)
\item make.plot - if to make a plot of the final inferred network (default: FALSE)
\item scale - the scaling parameter for color-coding (default: 1.5)
\item indepTest - the independence test wrapper function (default: mcX2Test). 
The following functions are available: \code{mcX2Test} (a wrapper around mc-x2-c (Monte Carlo X2 test)
with B=5000), \code{mcX2TestB50k} (a wrapper around mc-x2-c (Monte Carlo X2 test) test with B=50000), 
\code{mcMITest} (wrapper around mc-mi test from \code{bnlearn} with B=5000). 
The package \code{pcalg} additionally provide following tests:
\code{binCItest} for binary data (performs a G^2 test) and \code{gaussCItest} for continuous data 
(performs Fisher's Z transformation), \code{dicCItest} for discrete data (performs G^2 test). 
\item class.label - the label to show for target variable (default: "target")
\item use.colors - if to colour code the results (default: TRUE)
}}
}
\examples{# load in the data for fly mesoderm
data(mesoBin)

# increase alpha to 0.1, suppress progress output
calcDependence(mesoBin$VM, "ncpc", alpha=0.05)

# run ncpc* with mutual information with shrinkage and minimal numbers of 
# samples per conditioning set of 15
calcDependence(mesoBin$VM, "ncpc*", test.type="mi-sh", min.table.size=15)

# run PC algorithm using the G^2 test from pcalg package
calcDependence(mesoBin$VM, "pc", indepTest=pcalg::binCItest)

# run hill-climbing with BIC penalization and plot the resulting Bayesian Network
# NOTE: plotting requires the Rgraphviz package
if(require("Rgraphviz"))
calcDependence(mesoBin$VM, "hc-bic", make.plot=TRUE)

# continuous data example
data(mesoCont)

# run ncpc with linear correlation test and with maximal conditioning set of 3
res <- calcDependence(mesoCont$VM, "ncpc", max.set.size=3, test.type="cor")
# plot the resulting ddgraph with colours
if(require("Rgraphviz"))
plot(res$obj, col=TRUE)}

