\name{color.legend.DDGraph}
\alias{color.legend.DDGraph}
\title{Plot color coding legend}
\usage{color.legend.DDGraph(xl, yb, xr, yt, legend, rect.col, cex=1, align="lt", gradient="x",
    title="", ...)
}
\description{This function is a slightly modified version of function \code{color.legend()} function from plotrix package. 
It plots a color legend at the given coordinates. This version extends the original plotrix function with
additional label and ability to plot into margins.}
\arguments{\item{xl}{lower left corner x coordinate}
\item{yb}{lower left corner y coordinate}
\item{xr}{upper right corner x coordinate}
\item{yt}{upper right corner y coordinate}
\item{legend}{the text to be plotted below the color coding rectangle}
\item{rect.col}{the color that will fill the rectangle}
\item{cex}{character expansion factor for the labels}
\item{align}{how to align the labels relative to the color rectangle}
\item{gradient}{whether to have a horizontal (x) or vertical (y) color gradient}
\item{title}{the title to be printed above the color coding rectangle}
\item{...}{the additional arguments passed to \code{text()}}
}

