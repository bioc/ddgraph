\name{mapEnrichmentToColors}
\alias{mapEnrichmentToColors}
\title{Map enrichment values to colors...}
\usage{mapEnrichmentToColors(obj, palette, class.col, scale="auto")
}
\description{Map enrichment values to colors}
\details{The enrichment of every variable is calculated during construction of DDGraph objects (in \code{ncpc()}). 
Use this information to color code the node in the graph. By default the Orange-Red is used and shown the strength
of enrichment and depletion. No difference is made for enriched/depleted variables.}
\value{the p values color-coded by convertPvalueToColorIndex() function}
\arguments{\item{obj}{an object of type DDGraph}
\item{palette}{the color palette to use (by default Orange-Red)}
\item{class.col}{the color to use for class labels, if applicable (by default light green)}
\item{scale}{by how much to scale the -log10(p.value) when color coding: either a number of "auto" for automatic}
}

