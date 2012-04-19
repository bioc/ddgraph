#
# S4 generics
#

# FurlongDataSet
setGeneric("toDDDataSet", function(obj,...) standardGeneric("toDDDataSet"))
setGeneric("signalMatrix", function(object) standardGeneric("signalMatrix"))
setGeneric("classLabels", function(object) standardGeneric("classLabels"))

# DDDataSet
setGeneric("datasetName", function(obj, ...) standardGeneric("datasetName"))
setGeneric("rawData", function(obj, ...) standardGeneric("rawData"))
setGeneric("dataType", function(obj, ...) standardGeneric("dataType"))
setGeneric("variableNames", function(obj, ...) standardGeneric("variableNames"))
setGeneric("ciTest", function(obj,...) standardGeneric("ciTest"))

# DDGraph
#if(!isGeneric("plot"))
#	setGeneric("plot", function(x, y, ...) standardGeneric("plot"))
