\name{mapEnrichmentToColorsDual}
\alias{mapEnrichmentToColorsDual}
\title{Map enrichment values into two different palettes for enriched/depleted variables...}
\usage{mapEnrichmentToColorsDual(obj, palette.pos, palette.neg, class.col, scale="auto")
}
\description{Map enrichment values into two different palettes for enriched/depleted variables}
\value{the p values color-coded by convertPvalueToColorIndex() function}
\arguments{\item{obj}{an object of type DDGraph}
\item{palette.pos}{the palette to use for enrichment (by default Orange-Red)}
\item{palette.neg}{the palette to use for depletion (by default Purple-Blue)}
\item{class.col}{the colour to use for class labels, if applicable (by default light green)}
\item{scale}{by how much to scale the -log10(p.value) when color coding}
}
\examples{\dontrun{
data(mesoBin)
meso <- ncpc(mesoBin$Meso)
# use heat colours for both enrichment and depletion
mapEnrichmentToColorsDual(meso, palette.pos=heat.colors(10), palette.neg=heat.colors(10))
}}

