
setValidity("DDGraphEdge", function(object){
	if(object@type %in% c("directed", "undirected", "bidirectional", "dashed", "dotted"))
		TRUE
	else 
		paste("Invalid value for slot \"type\": \"", object@type, "\". Allowed values are: directed, undirected, bidirectional, dashed, dotted", sep="")
})

#' show method for DDGraphEdge
#' @param object the DDGraphEdge object
setMethod("show", signature=signature(object="DDGraphEdge"), function(object){
	edge.str = " "
	if( object@type == "directed" )
		edge.str = "->"
	if( object@type == "undirected" )
		edge.str = "-"
	if( object@type == "bidirectional" )
		edge.str = "<->"
	if( object@type == "dashed" )
		edge.str = "- -"
	if( object@type == "dotted" )
		edge.str = ".."
		
		
	cat("EDGE ", object@fromName, "[", object@fromInx, "] ", edge.str, " ", object@toName, "[", object@toInx, "] based on ", length(object@ciTests), " test(s)\n", sep="")
})
