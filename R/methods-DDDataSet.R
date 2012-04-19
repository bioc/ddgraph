#' Return the raw data frame with the variables, and the last column being "class"
#'
#' @title Raw data.frame with data
#' @param obj the DDDataSet object
#' @param ... unused
#' @return the raw dataframe that contains all the data
#' @export
setMethod("rawData", signature=signature(obj="DDDataSet"), function(obj,...) obj@data)

#' Return the data type ("binary" or "continuous")
#'
#' @title Return data type
#' @param obj the DDDataSet object
#' @param ... unused
#' @return the data type
#' @export
setMethod("dataType", signature=signature(obj="DDDataSet"), function(obj,...) obj@dataType)

#' Dataset name
#'
#' @param obj the DDDataSet object
#' @param ... unused
#' @return the name of the dataset used in plotting etc
#' @export
setMethod("datasetName", signature=signature(obj="DDDataSet"), function(obj,...) obj@name)

#' Names of variables (including "class")
#'
#' @title Names of variables (+class)
#' @param x the DDDataSet object
#' @return the names of the variables
#' @export
setMethod("names", signature=signature(x="DDDataSet"), function(x) colnames(x@data))

#' Names of variables (without "class")
#'
#' @title Names of variables (-class)
#' @param obj the DDDataSet object
#' @param ... unused
#' @return only the names of the variables (i.e. without "class")
#' @export
setMethod("variableNames", signature=signature(obj="DDDataSet"), function(obj, ...) colnames(obj@data)[1:(ncol(obj@data)-1)] )

#' access a specific variable in the dataset by name
#' @aliases $,DDDataSet-method
#' @param x the DDDataSet object
#' @param name the variable name
#' @docType methods
#' @rdname operators-DDDataSet
setMethod("$", signature=signature(x="DDDataSet"), function(x, name) x@data[[name]])

#' access a specific variable in the dataset by name
#' @aliases [[,DDDataSet-method
#' @param x the DDDataSet object
#' @param i variable name
#' @param j unusued
#' @rdname operators-DDDataSet
setMethod("[[", signature=signature(x="DDDataSet"), function(x, i, j) x@data[[i]])

#' access a specific variable in the dataset by name
#' @aliases [,DDDataSet-method
#' @param x the DDDataSet object
#' @param i rows
#' @param j columns
#' @param drop unused
#' @param ... unused
#' @rdname operators-DDDataSet
setMethod("[", signature=signature(x="DDDataSet", i="ANY", j="ANY"), function(x, i, j, ..., drop = TRUE) x@data[i,j])

#' Construct new DDDataSet object
#'
#' Try to initialise with anything that can be converted to matrix and vectors.
#' @param .Object the DDDataSet object
#' @param data the data slot
#' @param name the name slot
#' @param ... unused
setMethod("initialize", signature=signature(.Object="DDDataSet"), function(.Object, ..., data=data.frame(), 
          name=paste("Empty name created at", date())) {
	.Object@data = data
	.Object@name = name
	
	# check the types of variables and then decide on the type
	
	if( all(dim(data) == c(0,0))){
		.Object@dataType = "binary"
	}
	else if( is.binary(data) ){
		.Object@dataType = "binary"
	} else {
		.Object@dataType = "continuous"
	}
  
	.Object
})

#' show method for DDDataSet
#' @param object the DDDataSet object
setMethod("show", signature=signature(object="DDDataSet"), function(object){
	cat("DDDataSet object:", object@name, "\n")
	cat("with", nrow(object@data), "data points of", ncol(object@data)-1, "variables ")
	cat("with", object@dataType, "values\n")
})


#' Construct an DDDataSet object
#' 
#' @param signal the matrix or data frame where rows are observations and columns variables
#' @param name the name of the dataset (to be used in plotting, etc)
#' @param classLabels the vector of class labels or target responses (aka target variable)
#' @param classLabelsCol the column which should be interpreted as class labels (either name or index)
#' @param removeZeroVar if to remove zero variance columns without producing an error (default: TRUE)
#' 
#' @return a new DDDataSet object
#' @export
#' @examples
#' # columns are features, rows observations
#' data <- matrix(rbinom(50, 1, 0.5), ncol=5)
#' # target class labels
#' labels <- c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1)
#' makeDDDataSet(data, name="example data", classLabels=labels)
makeDDDataSet = function(signal, name, classLabels=NULL, classLabelsCol=NULL, removeZeroVar=FALSE){
  if(is.null(classLabels)){
    # check if column has been specified
    if(is.null(classLabelsCol)){
      if( !any(colnames(signal) == "class" ))
        stop("Either provide parameters classLabels or classLabelsCol, or name one of your input columns \"class\"")
      
      # extract and remove the class labels
      classLabels = signal[,"class"]
      signal = signal[, -which(colnames(signal) == "class")]
      
    } else {
      # extract the class labels column
      classLabels = signal[,classLabelsCol]
      # delete the column from signal
      if(is.numeric(classLabelsCol))
        signal = signal[,-classLabelsCol]
      else
        signal = signal[,-which(colnames(signal)==classLabelsCol)]
    }
  }
      
  if( is.vector(signal) ){
  	# convert to matrix
  	signal = matrix(signal, ncol=1)
  }
  
  # check if all column are numeric
  if( !all(apply(signal, 2, is.numeric)) )
    stop("All columns need to be numeric")
  if(!is.numeric(classLabels))
    stop("Class labels need to be numeric")
    
	if( is.binary(signal) & !is.binary(classLabels) )
		stop("Error making DDDataSet: class labels need to be binary (0,1) if all of data is binary")
  
  column.variance = apply(signal, 2, var)
  
  if(any(column.variance == 0)){
    zero.var = which(column.variance == 0)
    if( !removeZeroVar ){
      if( length(zero.var) == 1)
        stop(paste("Column", colnames(signal)[zero.var], 
                 "has zero variance, use option removeZeroVar=TRUE to automatically remove it"))
      else
        stop(paste("Columns:", paste(colnames(signal)[zero.var], collapse=","), 
                 "have zero variance, use option removeZeroVar=TRUE to automatically remove them"))
    } else {
      signal = signal[,-zero.var]
    }
  }
  
  	# check for duplicate column names
  	if(length(unique(colnames(signal))) != length(colnames(signal))){
  		stop("Column names of input parameter \"signal\" need to be unique.")
  	}
  
	# check if some of columns are duplicates
	duplicate = matrix(FALSE, nrow=ncol(signal), ncol=ncol(signal))
	for(i in 1:ncol(signal)){
		for(j in 1:ncol(signal)){
			if(j > i){
				duplicate[i,j] = identical(signal[,i], signal[,j])
			}
		}
	}
	
	if(any(duplicate)){
		dup.inx = which(duplicate, arr.ind=TRUE)
		if(!is.null(colnames(signal)))
			dup = apply(dup.inx, 1, function(x) paste(colnames(signal)[x[1]], "<->", colnames(signal)[x[2]]))
		else
			dup = apply(dup.inx, 1, function(x) paste(x[1], "<->", x[2]))
		dup = paste(dup, collapse=", ")
		stop("Columns with identical values are not allowed, please remove duplicates! Following pairs of columns have identical values: ", dup)
	}
		
	d = cbind(signal, "class"=classLabels)
	new("DDDataSet", data=as.data.frame(d), name=name)
}




