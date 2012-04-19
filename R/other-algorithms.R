#
# Front-end for Bayesian Network inference algorithms from bnlearn and pcalg packages
#


#' A custom plotting function for the BNlearn graphs
#'
#' @param d an object of type DDDataSet
#' @param bnlearn.function.name the bnlearn reconstruction algorithm to use (default: hc)
#' @param alpha the alpha value of conditional independence tests (if applicable)
#' @param test the type of conditional independence test (if applicable)
#' @param make.plot if to make a plot or just return the network (default: FALSE)
#' @param blacklist a data frame with two columns (optionally labeled "from" and "to"), 
#'        containing a set of arcs not to be included in the graph.
#' @param B the number of bootstrap runs of permutations (for iamb and such algorithms)
#' @param restart the number of random restarts for score-based algorithms
#' @param scale the color scaling
#' @param class.label the label to use for the class variable
#' @param use.colors if to color code the results
#' @param score the scoring penalization metric to use (when applicable)
#' @return an object of class "bn" representing the inferred network
#' @export
#' @examples
#' data(mesoBin)
#' # use hill-climbing to make the causal network and plot with enrichment colours
#' plotBNLearn(mesoBin$Meso, make.plot=TRUE)
plotBNLearn = function(d, bnlearn.function.name="hc", alpha=0.05, test="mc-mi", 
                       make.plot=FALSE, blacklist=NULL, B=NULL, restart=0, scale=1.5,
                       class.label="target", use.colors=TRUE, score="bic"){
	# for discrete data bnlearn requires a factor representation
	if(dataType(d) == "binary")
		dd = convertToFactor(rawData(d))
	else
		dd = rawData(d)
 
	colnames(dd)[colnames(dd)=="class"] = class.label
	
	# use one of the bnlearn functions
	bnlearn.function = eval(parse(text=bnlearn.function.name))
	if(bnlearn.function.name == "hc" | bnlearn.function.name == "mmhc")
		bn <- bnlearn.function(dd, blacklist=blacklist, restart=restart, score=score)
	else
		bn <- bnlearn.function(dd, test=test, alpha=alpha, blacklist=blacklist, B=B)
		
	if(make.plot){
	    if(!require("Rgraphviz")){
	    	stop("This function requires RGraphviz")
	    }
		# setup the highlight parameters
		highlight = list()
		# find the color mapping
		if(use.colors){
			eg = ncpc(d, max.set.size=0)
			ret.col = mapEnrichmentToColorsDual(eg, scale=scale)
			col = ret.col$col
			names(col)[names(col)=="class"] = class.label
			highlight[["nodes"]] = names(col) # we modify all nodes
			highlight[["fill"]] = col # set the fill colors
			highlight[["col"]] = rep("black", length(col)) # edge is still black
		} else {
			highlight = NULL
		}
	
		# plot
		graphviz.plot(bn, shape="ellipse", highlight=highlight, main=paste(toupper(bnlearn.function.name), 
                      ": ", datasetName(d), sep=""))
	}
	bn
}




#' Custom plotting function of PC algorithm to have nice highlighting
#'
#' @title Custom plotting for pcalgo
#' @param x an object of one of the pcalg classes
#' @param main the main title
#' @param labels the labels of each of the nodes (doesn't need to be a named vector)
#' @param colors the colors we want to assign to each node (doesn't need to be a named vector) 
#' @param ... additional parameters to pass to \code{layoutGraph()}
#' @export
customPlotPCAlgo = function (x, main = NULL, labels = NULL, colors = NULL, ...){
    if (is.null(main)) 
        main <- deparse(x@call)
    attrs <- list()
    nodeAttrs <- list()
    # if there are labels show the nodes as ellipses
    if (!is.null(labels)) {
        attrs$node <- list(shape = "ellipse", fixedsize = FALSE)
        names(labels) <- nodes(x@graph)
        nodeAttrs$label <- labels        
    }
    # layout the graph so we can change node colors, etc..
    if(!require("Rgraphviz")){
    	stop("This function requires RGraphviz")
    }
    layout = layoutGraph(x@graph, attrs = attrs, nodeAttrs = nodeAttrs, ...)
    
    if(!is.null(colors)){
    	names(colors) = nodes(x@graph)
	    nodeRenderInfo(layout)[["fill"]] = colors
	}
	
	# render the graph with additional parameters (like the main title)
    renderGraph(layout, graph.pars=list(graph=list(main=main)))

}


#' Infer a network using PC algorithm and plot it. 
#'
#' @title Plot the network inferred by the PC algorithm
#' @param d DDDataSet object
#' @param name the name to show (defaults to dataset name)
#' @param alpha the alpha value cut-off for the conditional independence tests
#' @param verbose if to show progress
#' @param directed if TRUE applies PC algorithm, if FALSE applies PC-skeleton
#' @param make.plot if to output the plot into the active device
#' @param scale the scaling parameter for color-coding
#' @param indepTest the independence test wrapper function (as needed by package \code{pcalg})
#' @param class.label the label to show for class variable
#' @param use.colors if to color code the results
#' @export
#' @examples
#' data(mesoBin)
#' # use PC algorithm to construct a causal network and colour it according to enrichment/depletion
#' plotPCalg(mesoBin$Meso, alpha=0.05, directed=TRUE, make.plot=TRUE)
plotPCalg = function(d, name=NULL, alpha=0.05, verbose=FALSE, directed=TRUE, make.plot=FALSE, 
                     scale=1.5, indepTest=mcX2Test, class.label="target", use.colors=TRUE){
	# number of variables and dataset name
	p = ncol(rawData(d))
	if(is.null(name))
	  name = datasetName(d)
	
	dd = rawData(d)
 
  colnames(dd)[colnames(dd)=="class"] = class.label
	# define independence test and sufficient statistics
	#indepTest <- mcX2Test # binCItest
	suffStat = list(dm = dd, adaptDF = FALSE)
	
	if(make.plot){
		# find out colouring of nodes	
		if(use.colors){
			eg = ncpc(d, max.set.size=0)
			colors = mapEnrichmentToColorsDual(eg, scale=scale)$col
			names(colors)[names(colors)=="class"] = class.label
		} else{
			colors = NULL
		}
	}
 
	# calculate the graph 
	if(directed){
		# use PC algorithm for directed graphs 
		fit <- pc(suffStat, indepTest, p, alpha, verbose = verbose)
		if( make.plot )
			customPlotPCAlgo(fit, main = paste("PC:",name), labels=names(dd), colors=colors)
	} else{
		# use only skeleton reconstruction
		fit <- skeleton(suffStat, indepTest, p, alpha, verbose = verbose)
		if( make.plot )
			customPlotPCAlgo(fit, main = paste("PC Skeleton:",name), labels=names(dd), colors=colors)
	}
	
	fit
}





#' Find the markov blanket for the PC algorithm output
#' 
#' @param pc output of PC algorithm from package pcAlgo, object of class "pcAlgo"
#' @param node the index of the node for which we are seeking the Markov Blanket
#' @return the inidicies of nodes that constitute the Markov Blanket
#' @export
pcalgMB = function(pc, node){
	p = as(pc@graph, "graphAM")
	adj = p@adjMat
	
	parents = which(adj[,node] == 1)
	children = setdiff(which(adj[node,] == 1), parents)
	spouses = c()
	for(i in children){
		spouses = union(spouses, which(adj[,i] == 1))
	}
	
	setdiff(union(union(parents, children), spouses), node)
}

#' Find the neighbourhood for the PC algorithm output
#' 
#' @param pc output of PC algorithm from package pcAlgo, object of class "pcAlgo"
#' @param node the index of the node for which we are seeking the neighbourhood
#' @return the inidicies of nodes that constitute the (undirected) neighbourhood
#' @export
pcalgNBR = function(pc, node){
	p = as(pc@graph, "graphAM")
	adj = p@adjMat
	
	parents = which(adj[,node] == 1)
	children = which(adj[node,] == 1)
	
	setdiff(union(parents, children), node)
}
