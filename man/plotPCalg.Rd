\name{plotPCalg}
\alias{plotPCalg}
\title{Plot the network inferred by the PC algorithm}
\usage{plotPCalg(d, name, alpha=0.05, verbose=FALSE, directed=TRUE, make.plot=FALSE,
    scale=1.5, indepTest=mcX2Test, class.label="target",
    use.colors=TRUE)
}
\description{Infer a network using PC algorithm and plot it.}
\arguments{\item{d}{DDDataSet object}
\item{name}{the name to show (defaults to dataset name)}
\item{alpha}{the alpha value cut-off for the conditional independence tests}
\item{verbose}{if to show progress}
\item{directed}{if TRUE applies PC algorithm, if FALSE applies PC-skeleton}
\item{make.plot}{if to output the plot into the active device}
\item{scale}{the scaling parameter for color-coding}
\item{indepTest}{the independence test wrapper function (as needed by package \code{pcalg})}
\item{class.label}{the label to show for class variable}
\item{use.colors}{if to color code the results}
}
\examples{data(mesoBin)
# use PC algorithm to construct a causal network and colour it according to enrichment/depletion
plotPCalg(mesoBin$Meso, alpha=0.05, directed=TRUE, make.plot=TRUE)}

