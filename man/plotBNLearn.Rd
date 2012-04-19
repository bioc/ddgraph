\name{plotBNLearn}
\alias{plotBNLearn}
\title{A custom plotting function for the BNlearn graphs...}
\usage{plotBNLearn(d, bnlearn.function.name="hc", alpha=0.05, test="mc-mi",
    make.plot=FALSE, blacklist, B, restart=0, scale=1.5,
    class.label="target", use.colors=TRUE, score="bic")
}
\description{A custom plotting function for the BNlearn graphs}
\value{an object of class "bn" representing the inferred network}
\arguments{\item{d}{an object of type DDDataSet}
\item{bnlearn.function.name}{the bnlearn reconstruction algorithm to use (default: hc)}
\item{alpha}{the alpha value of conditional independence tests (if applicable)}
\item{test}{the type of conditional independence test (if applicable)}
\item{make.plot}{if to make a plot or just return the network (default: FALSE)}
\item{blacklist}{a data frame with two columns (optionally labeled "from" and "to"), 
containing a set of arcs not to be included in the graph.}
\item{B}{the number of bootstrap runs of permutations (for iamb and such algorithms)}
\item{restart}{the number of random restarts for score-based algorithms}
\item{scale}{the color scaling}
\item{class.label}{the label to use for the class variable}
\item{use.colors}{if to color code the results}
\item{score}{the scoring penalization metric to use (when applicable)}
}
\examples{data(mesoBin)
# use hill-climbing to make the causal network and plot with enrichment colours
plotBNLearn(mesoBin$Meso, make.plot=TRUE)}

