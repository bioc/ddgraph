\name{plot,DDGraph,missing-method}
\alias{plot,DDGraph,missing-method}
\alias{plot}
\title{Plot DDGraphs using RGraphviz...}
\usage{\S4method{plot}{DDGraph,missing}(x, y, ..., col, legend=FALSE, only.legend=FALSE, plot.class=TRUE,
    class.label=datasetName(x@dataset), ci.symbol="open",
    plot.pvals=TRUE, pvals.format=function(x) sprintf("\%.2f", x),
    pvals.fontsize=12, main)
}
\description{Plot DDGraphs using RGraphviz}
\arguments{\item{x}{DDGraph object}
\item{y}{unusued}
\item{col}{specifies the colors to be used to color nodes. Can be any of the following:
\itemize{
\item named vector of colors
\item logical value (TRUE = nodes colored in default 0.1 to 1e-3 range, FALSE = no node coloring)
\item list of parameters to pass to mapEnrichmentToColorsDual(), 
valid pameters are: "palette.pos", "palette.neg", "class.col", "scale", "max.color.index"
}}
\item{legend}{if to plot the color legend}
\item{only.legend}{if to plot only the legend}
\item{plot.class}{if to plot class labels node}
\item{class.label}{if plot.class=TRUE the label of the class node}
\item{plot.pvals}{if to plot p values on top of edges}
\item{ci.symbol}{the RGraphviz arrowtail/head symbol name for conditional independence tests}
\item{pvals.format}{a function to format the p values to be displayed on directed edges}
\item{pvals.fontsize}{the size of the font for p values}
\item{main}{main title}
\item{...}{other parameters passed to layoutGraph()}
}
\examples{\dontrun{
# load data
data(mesoBin)
# make DDGraph
g <- ncpc(mesoBin$Meso)

# default plot
plot(g)

# use colours
plot(g, col=TRUE)
}}

