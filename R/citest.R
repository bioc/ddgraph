#
# Conditional independence test functions
#


#' Do conditional independence test on DDDataSet
#'
#' This function does a conditional independence var1 indep var2 | cond. The following
#' test types are available (implemented by package \code{bnlearn}).
#'
#' For binary data:
#' \itemize{
#'   \item "fisher" - Fisher's exact test (only for unconditional independence)
#'   \item "mi" - Mutual Information (discrete)
#'   \item "mi-sh" - Mutual Information (discrete, shrinkage)
#'   \item "mc-mi" - Mutual Information (discrete, Monte Carlo)
#'   \item "aict" - AIC-like Test
#'   \item "x2" - Pearson's X^2
#'   \item "mc-x2" - Pearson's X^2 (Monte Carlo)
#'   \item "mc-x2-c" - Pearson's X^2 (Monte Carlo) the corrected version
#'   \item "g2" - G^2 test (requires pcalg package)
#' }
#' For continuous data:
#' \itemize{
#'   \item "mi-g" - Mutual Information (Gaussian)
#'   \item "mi-g-sh" - Mutual Information (Gaussian, shrinkage)
#'   \item "mc-mi-g" - Mutual Information (Gaussian, Monte Carlo)
#'   \item "cor" - Pearson's Linear Correlation
#'   \item "mc-cor" - Pearson's Linear Correlation (Monte Carlo)
#'   \item "zf" - Fisher's Z Test
#'   \item "mc-zf" - Fisher's Z Test (Monte Carlo)
#' }
#'
#'
#' @param obj DDDataSet object on which (conditional) independence test needs to be done
#' @param var1 the name or index of the first variable to be tested
#' @param var2 the name or index of the second variable
#' @param cond the names or indexes of variables to condition on (defaults to NULL)
#' @param test.type the type of statistical test (defaults to mc-x2)
#' @param B the number of replicates for MC-based tests (default to NULL)
#' @param min.table.size the minimal number of samples in a contingency table per conditioning set 
#'        (makes sense only for discrete data)
#' @param ... unused
#' @return CITestResult object with the result of the test
#' @export
#' @examples
#' data(mesoBin)
#' # test if tin_4.6 is independent of class labels
#' ciTest(mesoBin$Meso, "Tin 4-6h", "class")
#' # test if tin_4.6 is independent of class conditioned on twi_2.4
#' ciTest(mesoBin$Meso, "Tin 4-6h", "class", "Twi 2-4h")
#' # repeat the test using G2 asymptotic distribution
#' ciTest(mesoBin$Meso, "Tin 4-6h", "class", "Twi 2-4h", test.type="g2")
setMethod("ciTest", signature=signature(obj="DDDataSet"), function(obj, var1, var2, cond=NULL, 
          test.type="mc-x2-c", B=NULL, min.table.size=NULL, ...){
	# these need to be single variables
	val.a = obj[[var1]]
	val.b = obj[[var2]]
	
	# check if variables exist
	if(is.null(val.a)){
		stop("Variable ", var1, " not found in ", datasetName(obj))
	}	
	if(is.null(val.b)){
		stop("Variable ", var2, " not found in ", datasetName(obj))
	}
	
	if( !is.null(cond) ){
		val.c = obj[,cond] # can be multiple columns
	} else {
		val.c = NULL
	}
	
	# find out the variable indexes and names
	if(is.numeric(var1)){
		var1.inx = var1
		var1.name = names(obj)[var1]
	} else {
		var1.inx = which(names(obj) == var1)
		var1.name = var1
	}
	
	if(is.numeric(var2)){
		var2.inx = var2
		var2.name = names(obj)[var2]
	} else {
		var2.inx = which(names(obj) == var2)
		var2.name = var2
	}
	
	if( is.null(cond) ){
		cond.inx = vector(mode="numeric")
		cond.name = vector(mode="character")
	} else if( is.numeric(cond) ){
		cond.inx = cond
		cond.name = names(obj)[cond]
	} else {
		cond.inx = which(names(obj) == cond)
		cond.name = cond
	}
	
	# check if it is a reliable test or not
	if(dataType(obj) == "binary" & !is.null(min.table.size) & length(cond) > 0){
		dd = as.data.frame(cbind(val.a, val.b, val.c))
		stopifnot(ncol(dd) == length(cond)+2)
		t = table(dd)
		
		sizes = as.vector(apply(t, 3:ncol(dd), sum))
		if(any(sizes < min.table.size)){
			return(new("CITestResult", targetInx=var1.inx, targetName=var1.name, sourceInx=var2.inx, sourceName=var2.name, 
			condSetInx=cond.inx, condSetName=cond.name, pValue=0, testType=test.type, reliable=FALSE))
		}
	}
	
	if( test.type == "fisher"){	
		# fisher's exact test
		if( dataType(obj) != "binary"){
			stop("Trying to do Fisher's test on continuous data!")		
		}		
		if( !is.null(cond) ){
			stop("Trying to do conditional independence test using Fisher's exact test (i.e. cond is not NULL)")
		}
		
		p.val = fisher.test(table(val.a, val.b))$p.value
	} else if( test.type == "mc-x2-c" ){
		if( dataType(obj) != "binary" )
			stop("mc-x2-c test only works on binary data!")
			
		if(length(cond) == 1)
			p.val = myX2c(val.a, val.b, obj[,cond], B=B)
		else
			p.val = myX2c(val.a, val.b, as.list(as.data.frame(obj[,cond])), B=B)
	} else if( test.type == "g2" ){
		if( dataType(obj) != "binary" )
			stop("g2 test only works on binary data!")
		
		if( is.null(cond) ){
			suffStat = list(dm = cbind(val.a, val.b), adaptDF = FALSE)
			p.val = binCItest(1, 2, NULL, suffStat)
		} else {
			suffStat = list(dm = cbind(val.a, val.b, val.c), adaptDF = FALSE)
			p.val = binCItest(1, 2, 3:ncol(suffStat$dm), suffStat)
		}
	} else {
		# only pass B if it is an mc-type of test
		if( length( grep("^mc", test.type) ) == 0 )
			B = NULL 
		# assume bnlearn test if other test names
		if( dataType(obj) == "binary" ){
			p.val = bnlearn::ci.test(convertToFactor(val.a), convertToFactor(val.b), convertToFactor(val.c), test=test.type, B=B)$p.value
		} else {
			p.val = bnlearn::ci.test(val.a, val.b, val.c, test=test.type)$p.value
		}
			
	}
	
	new("CITestResult", targetInx=var1.inx, targetName=var1.name, sourceInx=var2.inx, sourceName=var2.name, 
		condSetInx=cond.inx, condSetName=cond.name, pValue=p.val, testType=test.type, reliable=TRUE)
})

#' Implements the mc-x2 test in format needed for pcalg.
#' 
#' @title Wrapper around the bnlearn mc-x2 test
#' @param x the index of the first variable 
#' @param y the index of the second variable
#' @param S the conditioning set
#' @param suffStat the sufficient statistics to do the test, in this case a list of one element: 
#'        dm where the values matrix is stored
#'
#' @return p value of the test
#' @export
#' @examples 
#' suffStat <- list(dm = cbind("a"=c(0,1,0,0,1,0), "b"=c(1,0,0,0,1,0), "c"=c(0,0,0,1,1,1)))
#' # test if a is independent of b
#' mcX2Test(1, 2, NULL, suffStat)
#' # test if a is independent of b conditioned on c
#' mcX2Test(1, 2, 3, suffStat)
mcX2Test = function (x, y, S, suffStat) {
    dm <- suffStat$dm
    
    a = dm[,x]
    b = dm[,y]
    if( is.null(S) | length(S) == 0){
    	c = NULL
    } else {
    	if( length(S) == 1)
    		c = dm[,S]
    	else
    		c = as.list(as.data.frame(dm[,S]))
	}
    
    pval <- myX2c(a, b, c)
    pval
}

#' Implements the mc-mi test in format needed for pcalg.
#' 
#' @title Wrapper around the bnlearn mc-x2 test
#' @param x the index of the first variable 
#' @param y the index of the second variable
#' @param S the conditioning set
#' @param suffStat the sufficient statistics to do the test, in this case a list of one element: 
#'        dm where the values matrix is stored
#'
#' @return p value of the test
#' @export
#' @examples 
#' suffStat <- list(dm = cbind("a"=c(0,1,0,0,1,0), "b"=c(1,0,0,0,1,0), "c"=c(0,0,0,1,1,1)))
#' # test if a is independent of b
#' mcMITest(1, 2, NULL, suffStat)
#' # test if a is independent of b conditioned on c
#' mcMITest(1, 2, 3, suffStat)
mcMITest = function (x, y, S, suffStat) {
    dm <- suffStat$dm
    
    a = convertToFactor(dm[,x])
    b = convertToFactor(dm[,y])
    if( is.null(S) | length(S) == 0)
    	c = NULL
    else
	    c = convertToFactor(dm[,S])
    
    pval <- bnlearn::ci.test(a, b, c, test="mc-mi")$p.value
    pval
}

#' Version of \code{mcX2Test()} with 50000 Monte Carlo replicates. 
#' 
#' @title Wrapper around the bnlearn mc-x2 test (B=50k)
#' @param x the index of the first variable 
#' @param y the index of the second variable
#' @param S the conditioning set
#' @param suffStat the sufficient statistics to do the test, in this case a list of one element: 
#'        dm where the values matrix is stored
#'
#' @return p value of the test
#' @export
#' @examples 
#' suffStat <- list(dm = cbind("a"=c(0,1,0,0,1,0), "b"=c(1,0,0,0,1,0), "c"=c(0,0,0,1,1,1)))
#' # test if a is independent of b
#' mcX2TestB50k(1, 2, NULL, suffStat)
#' # test if a is independent of b conditioned on c
#' mcX2TestB50k(1, 2, 3, suffStat)
mcX2TestB50k = function (x, y, S, suffStat) {
    dm <- suffStat$dm
    
    a = dm[,x]
    b = dm[,y]
    if( is.null(S) | length(S) == 0){
    	c = NULL
    } else {
    	if( length(S) == 1)
    		c = dm[,S]
    	else
    		c = as.list(as.data.frame(dm[,S]))
	}
    
    pval <- myX2c(a, b, c, B=50000)
    pval
}

#' the inner loop for myX2c is implemented in C
#'
#' @param B the number of Monte Carlo replicates
#' @param numTable the number of conditional tables
#' @param rowSums the matrix or row sums for each conditional table (numTables x 4)
#' @param colSums the matrix or column sums for each conditional table (numTables x 4)
#' @return The values of chi-square statistics from random runs
mcX2CLoop = function(B, numTable, rowSums, colSums) {
	ret = .Call("mcX2CLoopC", B, numTable, rowSums, colSums, PACKAGE="ddgraph")
	#ret = .Call("mcX2CLoopC", B, numTable, rowSums, colSums)
	return(ret)
}

#' Get the value of chi-square statistics
#'
#' @param x is the contingency table
#' @param correct if to do the Yates correction
#' @return chisq statistics
chisq.val = function(x, correct=FALSE){
        n <- sum(x)
        
        nr <- nrow(x)
        nc <- ncol(x)
	    sr <- rowSums(x)
        sc <- colSums(x)
        
        if(any(sr == 0 | sc == 0) )
        	return(0)
        
        E <- outer(sr, sc, "*")/n

        if (correct && nrow(x) == 2 && ncol(x) == 2) {
                YATES <- 0.5
        } else {
                YATES <- 0
        }
        
        stats <- sum((abs(x - E) - YATES)^2/E)
        stats
}

#' The Monte-Carlo chi-square test
#'
#' This is the reimplementation of Monte Carlo chi-square test to be sure
#' it works correctly. The Monte Carlo loop is implemented using \code{Rcpp} and uses
#' the R function \code{r2dtable()} to generate random contingency tables with
#' fixed marginals. 
#'
#' @param x the first variable (vector of values)
#' @param y the second variable (vector of values)
#' @param C the variables to condition on - either a vector, or a list of vectors
#' @param B the number of Monte Carlo runs (defaults to 5000 if given NULL)
#' @return the P-value of the test
myX2c = function(x, y, C=NULL, B=5000){
	if(is.null(B))
		B = 5000

	invars = list(x,y)
	if( !is.null(C) ){
		if(is.list(C))
			invars = c(invars, C)
		else
			invars[[3]] = C
	}
	t = table(invars)
	
	# vector representation of contingency table
	v = as.vector(t)
	num.table = length(v)/4	
	stat = 0
	for(i in 1:num.table){
		i1 = 4*(i-1)+1
		i2 = i1 + 3
		stat = stat + chisq.val(matrix(v[i1:i2], ncol=2, nrow=2))
	}
	# calculate row and column sums
	row.sums = matrix(0, nrow=num.table, ncol=2)
	col.sums = matrix(0, nrow=num.table, ncol=2)
	for(i in 1:num.table){
		i1 = 4*(i-1)+1
		row.sums[i,1] = v[i1]+v[i1+2]
		row.sums[i,2] = v[i1+1]+v[i1+3]
		col.sums[i,1] = v[i1]+v[i1+1]
		col.sums[i,2] = v[i1+2]+v[i1+3]
	}
	
	# now simulate
	res = mcX2CLoop(B, num.table, row.sums, col.sums)	
	# there might be numerical differences, so substract 1e-8
	p.val = sum(res >= (stat - 1e-8)) / B
	p.val
	
}


#' Provide a unique ID composing of target, source and conditioning set (all names)
#' 
#' @param citest a CITestResult object
#'
#' @return a character ID
CITestResultID = function(citest){
	if(class(citest) != "CITestResult")
		stop("Input argument citest need to be of class CITestResult")
		
	paste(citest@targetName, "vs", citest@sourceName, "|", paste(citest@condSetName, collapse=","))	
}

#' Return a string representation of a variable represented with this CITest
#' 
#' @param citest an object of class CITestResult
CITestResultVar = function(citest){
	if(class(citest) != "CITestResult")
		stop("Input argument citest need to be of class CITestResult")
	
	if(length(citest@condSetInx) == 0)
		paste(citest@targetName)
	else
		paste(citest@targetName, ":", paste(citest@condSetName, collapse=","), sep="")
}
