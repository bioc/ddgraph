#
# Functions to deal with finding d-seperation properties, like active paths and blocking paths
#


#' Version of blockingNodes() for DDGraphs
#'
#' @param obj DDGraph object
#' @param nodes the selected nodes
#'
#' @return same as blockingNodes(): a list with blocking nodes and minimal length to the target node: target node => blocked by => number of steps
#' @export
blockingVariables = function(obj, nodes){
	if(class(obj) != "DDGraph")
		stop("This function requires a DDGraph")
		
	final.calls = obj@stats$final.calls
	res = list()
	for(n in nodes){
		inx = which(final.calls$name == n)
		res[[n]] = list()
		ex.by = final.calls$explained.by[inx]
		if(ex.by != "" & final.calls$type[inx] != "joint")
			res[[n]][[ex.by]] = 1
	}
	
	res

}

#' Find all such nodes in neighbourhood of source node that are blocking at least one active path leading to another node
#' 
#' @param allPaths a list of active paths from a source node (as produced by activePaths())
#' @param nodes a vector of target nodes for which we are finding blocking nodes
#'
#' @return a list with blocking nodes and minimal length to the target node: target node => blocked by => number of steps
#' @export
blockingNodes = function(allPaths, nodes){
	res = list()
	for(node in nodes){
		blocking = list()
		for(i in 1:length(allPaths)){
			if(node %in% allPaths[[i]]){
				p = allPaths[[i]]
				begin.inx = 2
				end.inx = which(p == node)
				if(end.inx > begin.inx){
					# blocking node from neighbhood of source node
					bl = as.character(p[begin.inx])
					if(bl %in% names(blocking)){
						if(blocking[[bl]] > (end.inx-1))
							blocking[[bl]] = end.inx-1
					} else {
						blocking[[bl]] = end.inx-1
					}
					#blocking = union(blocking, p[begin.inx])
				}
			}		
		}
		res[[ as.character(node) ]] = blocking
	}
	res
	
}

#' Find all active paths in a (partially) directed graph
#'
#' @param graph the graph either in one of the package \code{graph} classes, or of class \code{bn} or \code{pcAlgo}
#' @param node the source node of the path (index not name)
#' @param nodeNames optionally specify node names which can be used to return those instead of indicies
#'
#' @return a list of active paths with node as its source
#' @export
activePaths = function(graph, node, nodeNames=NULL){
	
	if(class(graph) == "bn"){
		adj = amat(graph)
	} else if(class(graph) == "pcAlgo"){
		p = as(graph@graph, "graphAM")
		adj = p@adjMat
	} else if(class(graph) %in% c("graphAM", "graphBAM", "graphNEL")){
		p = as(graph, "graphAM")
		adj = p@adjMat
	} else if(class(graph) == "matrix" && is.binary(graph) && nrow(graph) == ncol(graph)){
		adj = graph
	} else {
		stop("Unsupported input class of graph: ", class(graph))
	}
	
	if(!is.numeric(node)){
		stop("The node needs to be a numeric value indicating the index of the node")
	}
	
	# recursively find all active paths 
	# 
	# adj - adjacency matrix, if adj[i,j]=1 then there is an edge i~j
	# trail - vector of visited nodes
	# last - last edge followed: in (into the current node), out (out of current node) or none (first node)
	# allPaths - the list with all finished paths
	recursivePath = function(adj, trail, last, allPaths){
		i = trail[length(trail)]
		# if we extended the trail by another step
		stepMade = FALSE
		for(j in 1:ncol(adj)){
			if( !(j %in% trail) ){
				# A -> C -> B and A <- C -> B
				if(adj[i,j] == 0 & adj[j,i] == 0){
					next
				} else if(adj[i,j] == 1 & adj[j,i] == 0){
					allPaths = recursivePath(adj, c(trail, j), "in", allPaths)
					stepMade = TRUE
				} else if(adj[i,j] == 0 & adj[j,i] == 1 & last != "in"){
					allPaths = recursivePath(adj, c(trail, j), "out", allPaths)
					stepMade = TRUE
				} else if(adj[i,j] == 1 & adj[j,i] == 1){
					allPaths = recursivePath(adj, c(trail, j), "out", allPaths)
					stepMade = TRUE
				}
			}
		}
		if(!stepMade)
			allPaths[[length(allPaths)+1]] = trail
			
		allPaths	

	}	
		
	allPaths = recursivePath(adj, node, "none", list())
	
	# convert node indicies to names
	if(!is.null(nodeNames)){
		for(i in 1:length(allPaths)){
			allPaths[[i]] = nodeNames[ allPaths[[i]] ]
		}
	}
	
	allPaths
}


