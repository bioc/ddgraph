#
# Furlong dataset related function - mainly to load in and binarize the dataset
#


#' Get the names of variables (column names of signal matrix)
#' @param x FurlongDataSet object
#' @title Names of variables
#' @export
setMethod("names", signature=signature(x="FurlongDataSet"), function(x) colnames(x@signalMatrix))

#' Retrieve the matrix with raw signal values 
#' @param object FurlongDataSet object
#' @title Raw values
#' @export
setMethod("signalMatrix", signature=signature(object="FurlongDataSet"), function(object) object@signalMatrix)

#' Retrieve the vector of class labels (as factors)
#' @param object FurlongDataSet object
#' @title Class labels
#' @export
setMethod("classLabels", signature=signature(object="FurlongDataSet"), function(object) object@targetClasses)

#' Make the DDDataSet objects by selecting different tissues
#' 
#' @title DDDataSet object from FurlongDataSet
#' @param obj the FurlongDataSet object
#' @param tissues tissue names for which DDDataSet objects should be generated (default to all available tissues)
#' @param convertToBinary if to convert the signal into binary values
#' @param prettyNames if to make the names pretty, e.g. twi_2.4 -> Twi 2-4h
#' @param ... unused
#'
#' @return either single DDDataSet object, or a list of them (depending on number of selected tissues)
#' @export
#' @examples
#' # load binarized data with prettified names
#' all.data <- toDDDataSet(readFurlongData(), prettyNames=TRUE)
#' # load continuous data with original names
#' all.data <- toDDDataSet(readFurlongData(), convertToBinary=FALSE)
setMethod("toDDDataSet", signature=signature(obj="FurlongDataSet"), function(obj, tissues=c(), 
          convertToBinary=TRUE, prettyNames=FALSE, ...){

	# if tissues list is empty, take them all
	if( length(tissues) == 0 ){
		tissues = levels(classLabels(obj))
	}

	ret = list()
	# loop over all tisuses
	for( tissue in tissues ){	
		signal = signalMatrix(obj)
		labels = as.numeric(classLabels(obj) == tissue)
		# if needed convert to binary signal
		if(convertToBinary)
			signal[signal>0] = 1
		else
			labels[labels==0] = -1
			
		if(prettyNames){
			# prettify the names a bit
			n = colnames(signal)
			for(i in 1:length(n)){
				m = chartr("_", " ", n[i])
				m = chartr(".", "-", m)
				m = paste(toupper(substr(m, 1, 1)), substr(m, 2, nchar(m)), "h", sep="")
				n[i] = m
			}
			colnames(signal) = n
		}
			
		# make the DDDataSet object
		ret[[tissue]] = makeDDDataSet(signal, chartr("_", "&", tissue), labels)
	}

	# return either one object, or a names list	
	if( length(ret) == 1){
		return( ret[[1]] )
	} else{
		return( ret )
	}
		
})


#' Read the Furlong data into a FurlongDataSet object.
#'
#' Read the Furlong Dataset form the Supplementary Table 8 file provided with the package. An alternative
#' filename can be specified as well. 
#'
#' @title Read the Furlong Dataset 
#' @param infile the filename to load from, default to supplementary_table_8_training_set.txt in extdata/ of package
#'
#' @return an object of type FurlongDataSet witht the loaded data
#' @export
#' @examples
#' # read the furlong dataset that is provided with the package
#' readFurlongData()
readFurlongData = function(infile=NULL){
	if(is.null(infile))
		infile = dir(system.file(package="ddgraph",dir="extdata"), 
                     full.names=T, pattern=glob2rx("supplementary_table_8_training_set.txt"))

	tab8 = read.delim(infile, comment.char="#" )

	# extract all of signal
	signal.cols = which(names(tab8)=="tin_2.4") : ncol(tab8)
	signal = tab8[,signal.cols]

	# extract enhancer class
	class.cols = which(names(tab8)=="Meso") : which(names(tab8)=="VM_SM")
	tissues = names(tab8)[class.cols]
	enh.class = tab8[,class.cols]

	# convert class to single vector
	all.class = factor(colSums(apply(enh.class, 1, function(x) x*(1:length(x)))), labels=c("neg", names(enh.class)))
	
	# return the corresponding object
	new("FurlongDataSet", signalMatrix=as.matrix(signal), targetClasses=all.class)
}



