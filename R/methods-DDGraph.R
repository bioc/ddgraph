#' Construct new DDGraph object
#'
#' Properly initialize the object
#' @param .Object DDGraph object
#' @param direct direct variable indexes
#' @param indirect indirect variable indexes
#' @param joint joint variable indexes
#' @param conditional conditional variable indexes
#' @param conditionalJoint conditionally joint variable indexes
#' @param edges edges list
#' @param dataset DDDataSet object
#' @param params parameters used to make this object
#' @param stats the statistics used to make this object
#' @param ... unused
setMethod("initialize",signature=signature(.Object="DDGraph"), function(.Object, ..., 
          direct=vector(mode="numeric"), indirect=vector(mode="numeric"), 
		  joint=vector(mode="numeric"), conditional=vector(mode="numeric"), 
		  conditionalJoint=vector(mode="numeric"), edges=list(), 
		  dataset=new("DDDataSet"), params=list(), stats=list()) {
	if(is.null(direct))
		direct = vector(mode="numeric")
		
	if(is.null(indirect))
		indirect = vector(mode="numeric")
		
	if(is.null(joint))
		joint = vector(mode="numeric")
		
	if(is.null(conditional))
		conditional = vector(mode="numeric")
		
	if(is.null(conditionalJoint))
		conditionalJoint = vector(mode="numeric")
	
	.Object@direct = direct
	.Object@indirect = indirect
	.Object@joint = joint
	.Object@conditional = conditional
	.Object@conditionalJoint = conditionalJoint
	.Object@edges = edges
	.Object@dataset = dataset
	.Object@params = params
	.Object@stats = stats
  
	.Object
})

#' Names of different properties that can be accessed with $ operator
#'
#' @title Names of properties
#' @param x the DDGataSet object
#' @return the names of the variables
#' @export
setMethod("names", signature=signature(x="DDGraph"), function(x){
	 c("params", "final.calls", "stats", "direct", "joint", "indirect", "conditional", "conditionalJoint", 
	 "directAndJoint", "directAndJointAndConditional") 
})


#' access a property by name
#' @aliases $,DDGraph-method
#' @param x the DDGraph object
#' @param name the variable name
#' @docType methods
#' @rdname operators-DDGraph
setMethod("$", signature=signature(x="DDGraph"), function(x, name){
	if(name == "directAndJoint"){
		c(x@direct, x@joint)
	} else if(name == "directAndJointAndConditional"){
		unique(c(x@direct, x@joint, x@conditional))
	} else if(name == "final.calls"){
		x@stats$final.calls
	} else if(name %in% c("direct", "joint", "indirect", "conditional", "conditionalJoint", "params", "stats")){
		slot(x, name)
	} else {
		stop("Uknown element ", name)
	}
})

#' show method for DDGraph
#' @param object the DDGraph object
setMethod("show", signature=signature(object="DDGraph"), function(object){
	algorithm = ifelse(object@params$star, "ncpc*", "ncpc")
	cat("DDGraph produced with", algorithm, "algorithm\n")
	cat("Direct:", names(object@dataset)[object@direct], "\n")
	cat("Joint:", names(object@dataset)[object@joint], "\n")
	if(object@params$star){
		cat("Conditional:", names(object@dataset)[object@conditional], "\n")
		cat("Conditional joint:", names(object@dataset)[object@conditionalJoint], "\n")
	}
	cat("Indirect:", names(object@dataset)[object@indirect], "\n")
	cat("Using P-value alpha cutoff = ", object$params$alpha, 
		" with conditional independence test = \"", object$params$test.type, "\"\n", sep="")
})

