#
# Misc functions, mostly helper and convenience functions used throughout the code
#


#' Check if a vector, data frame or matrix contains only binary (0,1) values.
#'
#' @title Check if data structure has binary data in it 
#'
#' @param x the input vector, data.frame or matrix
#' @return boolean TRUE or FALSE
#' @export
#' @examples
#' # works on vectors, matrices and data frames
#' is.binary(0)
#' is.binary(c(1, 0, 0, 1, 0))
#' is.binary(matrix(c(1,0), nrow=2, ncol=2))
#' is.binary(data.frame("a"=c(1,0), "b"=c(0,1)))
#'
#' # returns FALSE if not binary
#' is.binary(c(1, 2, 3))
is.binary = function(x){
	if(is.data.frame(x)){
		v = unique(as.vector(as.matrix(x)))
	} else if(is.matrix(x) | is.vector(x)){
		v = unique(as.vector(x))
	} else{
		stop("Unknown input class ", class(x), " to is.binary()")
	}
	
	if( is.numeric(v) ){
		all(v %in% c(0,1))
	} else{ 
		FALSE
	}
	
}

#' Convert a matrix, dataframe or vector into a factor representation. Each column is going
#' to be separately converted into a factor.
#'
#' @title Convert data to factor representation
#'
#' @param x the input vector, data.frame or matrix
#' @export
#' @examples
#' # works on vectors, matrices and data frames
#' convertToFactor(0)
#' convertToFactor(c(1, 0, 0, 1, 0))
#' convertToFactor(matrix(c(1,0), nrow=2, ncol=2))
#' convertToFactor(data.frame("a"=c(1,0), "b"=c(0,1)))
convertToFactor = function(x){
	# convert to factors column-wise
	if(is.null(x)){
		NULL
	} else if(is.matrix(x)){
		oldnames = colnames(x)
		out = data.frame( apply(x, 2, as.factor) )
		if( !all(is.null(oldnames)) ){
			colnames(out) = oldnames
		}
		out
	} else if(is.data.frame(x)){
		oldnames = colnames(x)
		out = data.frame( apply(x, 2, as.factor) )
		colnames(out) = oldnames
		out
	} else{
		as.factor(x)
	}
}

