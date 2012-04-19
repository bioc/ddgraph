#' Access slots using the dollar notation
#' @aliases $,CITestResult-method
#' @param x the CITestResult object
#' @param name the slot name
#' @docType methods
#' @rdname operators-CITestResult
setMethod("$", signature=signature(x="CITestResult"), function(x, name) slot(x, name))

#' Access slots using the double square bracket notation
#' @aliases [[,CITestResult-method
#' @param x the CITestResult object
#' @param i the slot name
#' @param j unused
#' @param ... unused
#' @rdname operators-CITestResult
setMethod("[[", signature=signature(x="CITestResult", i="ANY", j="ANY"), function(x, i, j, ...) slot(x, i))

#' Names of slots that can be accessed with $ notation
#' @param x the CITestResult object
#' @export
setMethod("names", signature=signature(x="CITestResult"), function(x) slotNames(x))

#' show method for CITestResult
#' @param object the CITestResult object
setMethod("show", signature=signature(object="CITestResult"), function(object){
	cat("Conditional independence test results using", object@testType, "\n")
	cat("  ", object@targetName, "[", object@targetInx, "] vs ", object@sourceName, "[", object@sourceInx, "], cond=", 
		paste(object@condSetName, collapse=","), "[", paste(object@condSetInx, collapse=","),  
		"], reliable=", object@reliable, ", p.value=", object@pValue, "\n", sep="")
})
