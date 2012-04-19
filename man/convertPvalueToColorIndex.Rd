\name{convertPvalueToColorIndex}
\alias{convertPvalueToColorIndex}
\title{Convert P-values to color index...}
\usage{convertPvalueToColorIndex(p.vals, scale="auto", max.color.index, minimal.p.value=1e-04)
}
\description{Convert P-values to color index}
\details{Convert p values to a color index to color nodes in a graph. The P-values are fit into a range from 1 to
\code{max.color.index} by applying a scale. Before fitting, P-value are transformed by taking a log10, 
and a minimal P-value is needed to avoid -Inf results for very small P-values. Scale can either be a number
or "auto" in which case color coding is such that all P-values fit into the range.}
\value{a \code{list} with following elements: \code{col} - the color indexes, \code{zlim} - the
actual scale range (in log10) over the colors}
\arguments{\item{p.vals}{the P-values (after any multiple testing correction)}
\item{scale}{the color is calculated liked -log10( p.value) * scale, thus scale is used to scale the -log10 to 
the desired range. Either a number or "auto" for automatic}
\item{max.color.index}{the maximal color index to return}
\item{minimal.p.value}{the minimal P-value we accept (since from Monte Carlo we can get 0)}
}
\examples{# scale the P values into the log10 space of [1e-3,1] represented by max 6 colours
convertPvalueToColorIndex(c(0.01, 0.2, 0.3), scale="auto", max.color.index=6, minimal.p.value=1e-3)}

