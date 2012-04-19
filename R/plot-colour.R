#
# Functions to deal with color-coding the plotted DDGraphs
#

#' Convert P-values to color index
#'
#' Convert p values to a color index to color nodes in a graph. The P-values are fit into a range from 1 to
#' \code{max.color.index} by applying a scale. Before fitting, P-value are transformed by taking a log10, 
#' and a minimal P-value is needed to avoid -Inf results for very small P-values. Scale can either be a number
#' or "auto" in which case color coding is such that all P-values fit into the range. 
#'
#' @param p.vals the P-values (after any multiple testing correction)
#' @param scale the color is calculated liked -log10( p.value) * scale, thus scale is used to scale the -log10 to 
#'        the desired range. Either a number or "auto" for automatic
#' @param max.color.index the maximal color index to return 
#' @param minimal.p.value the minimal P-value we accept (since from Monte Carlo we can get 0)
#'
#' @return a \code{list} with following elements: \code{col} - the color indexes, \code{zlim} - the
#'         actual scale range (in log10) over the colors
#' @export
#' @examples
#' # scale the P values into the log10 space of [1e-3,1] represented by max 6 colours
#' convertPvalueToColorIndex(c(0.01, 0.2, 0.3), scale="auto", max.color.index=6, minimal.p.value=1e-3)
convertPvalueToColorIndex = function(p.vals, scale="auto", max.color.index, minimal.p.value=1e-4){
	if( is.null(scale) ) 
		scale = "auto"
	
	# put it to minimal value
	if( any(p.vals < minimal.p.value) )
		p.vals[p.vals < minimal.p.value] = minimal.p.value
				
	color.vals = -log10( p.vals )
	
	# automatic scaling, scale colors.vals into [0, max.color.index]
	if( !is.numeric(scale) ){
		scale = max.color.index / max(color.vals)	
	}
	color.vals = color.vals * scale
	
	if(any(color.vals > max.color.index))
		color.vals[color.vals>max.color.index] = max.color.index
		
	list(col=round(color.vals), zlim=c(0, max.color.index/scale))
}

#' Map enrichment values to colors
#'
#' The enrichment of every variable is calculated during construction of DDGraph objects (in \code{ncpc()}). 
#' Use this information to color code the node in the graph. By default the Orange-Red is used and shown the strength
#' of enrichment and depletion. No difference is made for enriched/depleted variables. 
#'
#' @param obj an object of type DDGraph
#' @param palette the color palette to use (by default Orange-Red)
#' @param class.col the color to use for class labels, if applicable (by default light green)
#' @param scale by how much to scale the -log10(p.value) when color coding: either a number of "auto" for automatic
#'
#' @return the p values color-coded by convertPvalueToColorIndex() function
#' @export
mapEnrichmentToColors = function(obj, palette=NULL,	class.col=NULL, scale="auto"){
	mapEnrichmentToColorsDual(obj, palette, palette, class.col, scale)
}

#' Map enrichment values into two different palettes for enriched/depleted variables
#'
#' @param obj an object of type DDGraph
#' @param palette.pos the palette to use for enrichment (by default Orange-Red)
#' @param palette.neg the palette to use for depletion (by default Purple-Blue)
#' @param class.col the colour to use for class labels, if applicable (by default light green)
#' @param scale by how much to scale the -log10(p.value) when color coding
#'
#' @return the p values color-coded by convertPvalueToColorIndex() function
#' @export
#' @examples
#' \dontrun{
#' data(mesoBin)
#' meso <- ncpc(mesoBin$Meso)
#' # use heat colours for both enrichment and depletion
#' mapEnrichmentToColorsDual(meso, palette.pos=heat.colors(10), palette.neg=heat.colors(10))
#' }
mapEnrichmentToColorsDual = function(obj, palette.pos=NULL, palette.neg=NULL, class.col=NULL, scale="auto") {
	
	# our target pallette, make sure we have some values
	if(is.null(palette.pos))
		palette.pos = c("#FFFFFF", brewer.pal(9, "OrRd"))[1:7]
	if(is.null(palette.neg))
		palette.neg = c("#FFFFFF", brewer.pal(9, "PuBu"))[1:7]
	if(is.null(class.col))
		class.col = brewer.pal(9, "Greens")[3]
		
	if( length(palette.pos) != length(palette.neg) )
		stop("The positive and negative palettes (palette.pos, palette.neg) need to be of same length")
	
	# number of variables
	p.vals = obj@stats$setE.pvals
	enrich = obj@stats$setE.log2FC > 0
	
	color.vals.ret = convertPvalueToColorIndex( p.vals, scale=scale, max.color.index=length(palette.pos)-1 )
	color.vals = color.vals.ret$col	
	
	# select the colours from the pallete and add "lightblue" for class labels
	out1 = c( palette.pos[color.vals+1], class.col)
	out2 = c( palette.neg[color.vals+1], class.col)
	# combine them
	out = out1
	out[ !enrich ] = out2 [ !enrich ]
	names(out) = names(obj@dataset)
	
	if( all(palette.pos == palette.neg) ){
		list(col=out, zlim=color.vals.ret$zlim, palette=palette.pos)
	} else {
		list(col=out, zlim=c(-color.vals.ret$zlim[2], color.vals.ret$zlim[2]), 
			palette=c(rev(palette.neg[2:length(palette.neg)]), palette.pos) )
	}
}

#' This function is a slightly modified version of function \code{color.legend()} function from plotrix package. 
#' It plots a color legend at the given coordinates. This version extends the original plotrix function with
#' additional label and ability to plot into margins. 
#'
#' @title Plot color coding legend
#'
#' @param xl lower left corner x coordinate
#' @param yb lower left corner y coordinate
#' @param xr upper right corner x coordinate
#' @param yt upper right corner y coordinate
#' @param legend the text to be plotted below the color coding rectangle
#' @param rect.col the color that will fill the rectangle
#' @param cex character expansion factor for the labels
#' @param align how to align the labels relative to the color rectangle
#' @param gradient whether to have a horizontal (x) or vertical (y) color gradient
#' @param title the title to be printed above the color coding rectangle
#' @param ... the additional arguments passed to \code{text()} 
color.legend.DDGraph = function (xl, yb, xr, yt, legend, rect.col, cex = 1, align = "lt", 
    gradient = "x", title="", ...) {
	mult = 0.6	
	
    oldcex <- par("cex")
    par(xpd = NA, cex = cex)
    gradient.rect(xl, yb, xr, yt, col = rect.col, nslices = length(rect.col), 
        gradient = gradient)
    if (gradient == "x") {
        xsqueeze <- (xr - xl)/(2 * length(rect.col))
        textx <- seq(xl + xsqueeze, xr - xsqueeze, length.out = length(legend))
        if (match(align, "rb", 0)) {
            texty <- yb - mult * strheight("O")
            textadj <- c(0.5, 1)
        }
        else {
            texty <- yt + mult * strheight("O")
            textadj <- c(0.5, 0)
        }
    }
    else {
        ysqueeze <- (yt - yb)/(2 * length(rect.col))
        texty <- seq(yb + ysqueeze, yt - ysqueeze, length.out = length(legend))
        if (match(align, "rb", 0)) {
            textx <- xr + mult * strwidth("O")
            textadj <- c(0, 0.5)
        }
        else {
            textx <- xl - mult * strwidth("O")
            textadj <- c(1, 0.5)
        }
    }
    text(textx, texty, labels = legend, adj = textadj, ...)
    text(mean(c(xl,xr)), yt + mult * strheight("O"), labels = title, adj=c(0.5, 0), ...)
    par(xpd = FALSE, cex = oldcex)
}

