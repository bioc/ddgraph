#
# plot for DDGraph
#


#' Plot DDGraphs using RGraphviz 
#'
#' @param x DDGraph object 
#' @param y unusued
#' @param col specifies the colors to be used to color nodes. Can be any of the following:
#'        \itemize{
#'            \item named vector of colors
#'            \item logical value (TRUE = nodes colored in default 0.1 to 1e-3 range, FALSE = no node coloring) - only available for binary datasets.
#'            \item list of parameters to pass to mapEnrichmentToColorsDual(), 
#'                  valid pameters are: "palette.pos", "palette.neg", "class.col", "scale", "max.color.index"
#'        }
#' @param legend if to plot the color legend
#' @param only.legend if to plot only the legend
#' @param plot.class if to plot class labels node
#' @param class.label if plot.class=TRUE the label of the class node
#' @param plot.pvals if to plot p values on top of edges
#' @param ci.symbol the RGraphviz arrowtail/head symbol name for conditional independence tests
#' @param pvals.format a function to format the p values to be displayed on directed edges
#' @param pvals.fontsize the size of the font for p values
#' @param main main title
#' @param ... other parameters passed to layoutGraph()
#'
#' @export
#' @examples
#' \dontrun{
#' # load data
#' data(mesoBin)
#' # make DDGraph
#' g <- ncpc(mesoBin$Meso)
#'
#' # default plot
#' plot(g)
#'
#' # use colours
#' plot(g, col=TRUE)
#' }
setMethod("plot", signature=signature(x="DDGraph", y="missing"), function(x, y, ..., col=NULL, 
          legend=FALSE, only.legend=FALSE, plot.class=TRUE, class.label=datasetName(x@dataset), 
          ci.symbol = "dot", plot.pvals=TRUE, pvals.format=function(x) sprintf("%.2f", x), 
          pvals.fontsize=12, main=NULL){

	# make sure we have Rgraphviz          
    if (!require("Rgraphviz"))
	    stop("this function requires Rgraphviz.")
	    
	# x is the object, y is ignored
	obj = x
	
	col.zlim = NULL
	col.palette = NULL
	
	if(!is.null(class.label)){
		if(length(grep("\\|", class.label))> 0)
			class.label = chartr("|", " ", class.label)
		
		# make sure the dataset is named different than a variable	
		if(class.label %in% variableNames(obj@dataset)){
			stop("The dataset is called like one of the variable, please use the class.label parameter to override the dataset name in the plot")
		}
	}
	
	# parse the col parameter
	if( is.logical(col) ){
		if(col && all(!is.na(obj@stats$setE.log2FC))){
			ret.col = mapEnrichmentToColorsDual(obj)
			col = ret.col$col
			col.zlim = ret.col$zlim
			col.palette = ret.col$palette
		} else {
			col = NULL
		}
	} else if( is.list(col) ){
		valid.params = c("palette.pos", "palette.neg", "class.col", "scale")
		if( !all(names(col) %in% valid.params) )
			stop("Invalid color parameter supplied, valid parameters are:", paste(valid.params, ", "))
	
		ret.col = mapEnrichmentToColorsDual(obj, palette.pos=col[["palette.pos"]], palette.neg=col[["palette.neg"]], 
			class.col=col[["class.col"]], scale=col[["scale"]] )
		col = ret.col$col
		col.zlim = ret.col$zlim
		col.palette = ret.col$palette
	}
	
	# only plot the legend and exit
	if( only.legend) {
	
		if( is.null(col.zlim) ){
			stop("Legends only work when color scheme is specified. Either with col=TRUE or col=list(...).")
		}
	
		par(mar=c(1,1,1,1))
		plot(NULL, xlim=c(0,1), ylim=c(0,1), type="n", axes=F, xlab="", ylab="")
		
		legend.text = c(paste("1e",round(col.zlim[1]),sep=""), "<- depleted | enriched ->", paste("1e-", round(col.zlim[2]), sep=""))
		color.legend.DDGraph(0.0, 0.45, 1, 0.55, 
			legend.text, col.palette, cex=0.8, align="rb", title="Marginal dependence P-values")
		
		return(invisible(NULL))
	}
	
	# set the correct label for class
	if(plot.class & !is.null(col)){
		names(col)[names(col)=="class"] = class.label
	}
	
	dd = rawData(obj@dataset)
	# only the input variables (without class labels)
	vars = dd[,1:(ncol(dd)-1)]
	edges = obj@edges
	
	sel = c()
	for(edge in edges){
		sel = union(sel, edge@fromInx)
		sel = union(sel, edge@toInx)
	}
	
	sel = union(sel, obj@stats$setE.selected)
	
	# not enriched nodes might occur in DDGraph* graphs
	not.enriched = c(obj@conditional, obj@conditionalJoint)
	
	# all node names
	if(plot.class)
		all.nodes = c(names(vars)[sel], class.label)
	else
		all.nodes = names(vars)[sel]
		
	g = new("graphNEL", nodes = all.nodes, edgemode = "undirected")
	
	# add all the edges present in our model
	for(edge in edges){
		g = addEdge(edge@fromName, edge@toName, g, 1)
	}
	
	# if we are plotting the class node then add edges to all direct/joint variables
	if(plot.class){
		for(x in union(names(obj@direct), names(obj@joint)))
			g = addEdge(x, class.label, g, 1)
	}
	
	# set the shape for direct
	if(length(obj@direct)>0){
		shape = rep("ellipse",length(obj@direct))
		names(shape) = names(vars)[obj@direct]
		nodeRenderInfo(g) <- list(shape = shape)
	}
	
	# shapes for indirect
	if(length(obj@indirect)>0){
		shape = rep("box", length(obj@indirect))
		names(shape) = names(vars)[obj@indirect]
		nodeRenderInfo(g) <- list(shape = shape)	
	
		# put in the border color
		cols = rep("white", length(obj@indirect))
		names(cols) = names(vars)[obj@indirect]
		nodeRenderInfo(g) <- list(col = cols)
	}
	
	# shapes for joint
	if(length(obj@joint) > 0){
		shape = rep("ellipse", length(obj@joint))
		names(shape) = names(vars)[obj@joint]
		nodeRenderInfo(g) <- list(shape = shape)
	}
	
	# shapes for not enriched	
	if(length(not.enriched) > 0){
		shape = rep("ellipse", length(not.enriched))
		names(shape) = names(vars)[not.enriched]
		nodeRenderInfo(g) <- list(shape = shape)
	}

	
	if(plot.class){
		nodeRenderInfo(g) <- list(shape = structure("circle", names=class.label))
	}
	
	# set the border	
	all.lty = rep("solid", length(all.nodes))
	names(all.lty) = all.nodes
	
	if(length(obj@joint) > 0)
		all.lty[all.nodes %in% names(vars)[obj@joint]] = "dashed"
	if(length(not.enriched) > 0)
		all.lty[all.nodes %in% names(vars)[not.enriched]] = "dotted"
		
	all.lwd = rep(2, length(all.nodes))
	names(all.lwd) = all.nodes
		
	nodeRenderInfo(g) <- list(lty = all.lty, lwd=all.lwd)
	
	# put p values as labels if neccessary
	labels = list()
	if(plot.pvals){
		# add all edges
		for(edge in edges){
			if(edge@type %in% c("directed", "bidirectional")){
				edge.name = paste(edge@fromName, edge@toName, sep="~")
				if(!(edge.name %in% edgeNames(g)))
					edge.name = paste(edge@toName, edge@fromName, sep="~")
			
				# construct the label from all ciTests
				lab = ""
				for(ci in edge@ciTests){
					if(lab == "")
						lab = paste(pvals.format(ci@pValue))
					else
						lab = paste(lab,"/", pvals.format(ci@pValue), sep="")
				}
				
				labels[[edge.name]] = lab
			}
		}
		
	}
	

	################ layout the graph #############################

	if(plot.pvals && length(labels)>0){
		gl = layoutGraph(g, nodeAttrs=list(shape="ellipse"), edgeAttrs = list(label = labels), ...)
	} else {	
		gl = layoutGraph(g, nodeAttrs=list(shape="ellipse"), ...)
	}
	
	if(!is.null(col)){    	
		# the fill colours
	    nodeRenderInfo(gl)[["fill"]] = col
	} else{
		if(length(obj@indirect) > 0){
			cols = rep("#F0F0F0", length(obj@indirect))
			names(cols) = names(vars)[obj@indirect]
			nodeRenderInfo(gl) <- list("fill" = cols)
		}
	}
	
	# plot edges
	
	#arrow.symbol <- function(x, ...) {
	#	points(x, pch=ci.pch, cex=ci.cex, ...)
	#}
	
	arrow.symbol = ci.symbol
	for(edge in edges){	
		edge.name = paste(edge@fromName, edge@toName, sep="~")
	
		if(edge@type == "directed"){
			if(edge.name %in% edgeNames(g)){
				param = list(arrow.symbol)
				names(param) = edge.name
				edgeRenderInfo(gl) = list(arrowhead = param)
			} else{
				edge.name = paste(edge@toName, edge@fromName, sep="~")
				param = list(arrow.symbol)
				names(param) = edge.name
				edgeRenderInfo(gl) = list(arrowtail = param)
			}
		} else if(edge@type == "bidirectional"){
			if(!(edge.name %in% edgeNames(g)))
				edge.name = paste(edge@toName, edge@fromName, sep="~")
			param = list(arrow.symbol)
			names(param) = edge.name
			param2 = "both"
			names(param2) = edge.name
			param3 = "gray"
			names(param3) = edge.name
			edgeRenderInfo(gl) = list(arrowhead = param, arrowtail=param, direction = param2,
				col=param3)
		} else if(edge@type == "dashed"){
			if(!(edge.name %in% edgeNames(g)))
				edge.name = paste(edge@toName, edge@fromName, sep="~")
			param = "dashed"
			names(param) = edge.name
			param3 = "black"
			names(param3) = edge.name
			edgeRenderInfo(gl) = list(lty = param, col=param3)
		} else if(edge@type == "dotted"){
			if(!(edge.name %in% edgeNames(g)))
				edge.name = paste(edge@toName, edge@fromName, sep="~")
			param = "dotted"
			names(param) = edge.name
			param3 = "black"
			names(param3) = edge.name
			edgeRenderInfo(gl) = list(lty = param, col=param3)
		} else if(edge@type == "undirected"){
			if(!(edge.name %in% edgeNames(g)))
				edge.name = paste(edge@toName, edge@fromName, sep="~")
			param3 = "black"
			names(param3) = edge.name
			edgeRenderInfo(gl) = list(col=param3)
		}
	}	
	
	# properly colour edges leading to class label as well
	if(plot.class){
		for(x in union(names(obj@direct), names(obj@joint))){
			edge.name = paste(x, class.label, sep="~")
			if(!(edge.name %in% edgeNames(g)))
				edge.name = paste(class.label, x, sep="~")
				
			if(x %in% names(obj@joint))
				param3 = "gray"
			else
				param3 = "black"
			names(param3) = edge.name
			edgeRenderInfo(gl) = list(col=param3)
		}
			
	}
	
	# some global plotting parameters 
	if(is.null(main))
		main.title = paste("DDGraph for", datasetName(obj@dataset))
	else
		main.title = main
	
	graph.par(list(graph = list(main=main.title)))
	
	if( !is.null(col.zlim) & legend ){
		par(oma=c(6,0,0,0))
	}
		
	# draw
	graph.par(list(edges=list(fontsize=pvals.fontsize)))
	try({renderGraph(gl)})
	
	# add the colour legend
	if( !is.null(col.zlim) & legend){
		plot.space = par()[["usr"]]
		x1 = plot.space[1]
		x2 = plot.space[2]
		y1 = plot.space[3]
		y2 = plot.space[4]
		
		# FIXME: need better way to do this... 
		legend.text = c(paste("1e",round(col.zlim[1]),sep=""), "<- depleted | enriched ->", paste("1e-", round(col.zlim[2]), sep=""))
		color.legend.DDGraph(x1+(x2-x1)*(0.25), y1+(y2-y1)*(-0.32), x1+(x2-x1)*0.75, y1+(y2-y1)*(-0.27), 
			legend.text, col.palette, cex=0.8, align="rb", title="Marginal dependence P-values")
	}	

})
