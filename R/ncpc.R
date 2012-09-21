#
# NCPC/NCPC* algorithm implementation
#


#' Make a Direct Dependence Graph using the NCPC algorithm
#'
#' Make a Direct Dependence Graph using a P-value and conditional independence tests. There are two version of the
#' algorithm: NCPC and NCPC*. NCPC finds the causal neighbourhood while the NCPC* infers the full Markov Blanket. 
#'
#' The full algorithm is given in (Stojnic et al, 2012).
#' 
#' @param obj DDDataSet object
#' @param alpha the alpha (P-value) cutoff for conditional independence tests 
#' @param p.value.adjust.method the multiple testing correction adjustment method 
#' @param test.type the type of conditional independence test (default: Monte Carlo x2 test "mc-x2-c" for binary data
#'        and partial correlation "cor" for continuous data) . See the documentation 
#'        for \code{\link{ciTest}} for other available conditional independence tests
#' @param max.set.size the maximal number of variables to condition on, if NULL
#'        estimated from number of positives in class labels (default: NULL)
#' @param mc.replicates the number of Monte-carlo replicates, if applicable (default: 5000)
#' @param report.file name of the file where a detailed report is to be printed, 
#'        reporting is suppressed if NULL (default: NULL)
#' @param verbose if to print out information about how the algorithm is progressing (default: TRUE)
#' @param star if to use the NCPC* algorithm (default: FALSE)
#' @param min.table.size the minimal number of samples in a contingency table per conditioning set 
#'        (makes sense only for discrete data)
#'
#' @return DDGraph object	
#' @references R. Stojnic et al (2012): "A Graphical Modelling Approach to the Dissection of Highly Correlated Transcription Factor Binding Site Profiles", in press, PloS Computational Biology.
#' @export
#' @examples
#' 
#' ### load binary data for Mesoderm
#' data(mesoBin)
#' # run the NCPC algorithm with alpha=0.05 (on discrete data)
#' ncpc(mesoBin$Meso, alpha=0.05, test.type="mc-x2-c")
#' # run the NCPC* algorithm with alpha=0.05 (on discrete data)
#' res <- ncpc(mesoBin$Meso, alpha=0.05, test.type="mc-x2-c", star=TRUE)
#'
#' # analysis of results:
#' class(res)
#' # although of class DDGraph, behaves much like a list
#' names(res)
#' # parameters used in obtaining results
#' res$params
#' # labels for each of the variables
#' res$final.calls
#' # direct variables
#' res$direct
#'
#' ### load continous data
#' data(mesoCont)
#' # run the NCPC algorith with alpha=0.05 (on continuous data)
#' ncpc(mesoCont$Meso, alpha=0.05, test.type="cor", max.set.size=1)
#' # run the NCPC* algorith with alpha=0.05 (on continuous data)
#' ncpc(mesoCont$Meso, alpha=0.05, test.type="cor", max.set.size=1, star=TRUE)
#' 
ncpc = function(obj, alpha=0.05, p.value.adjust.method="none", test.type=c("mc-x2-c", "cor"), max.set.size=NULL, mc.replicates=5000, report.file=NULL, verbose=FALSE, star=FALSE, min.table.size=10){
    
    if(class(obj) != "DDDataSet"){
    	stop("obj needs to be of type DDDataSet")
    }
    
    # different default 
    if( identical(test.type, c("mc-x2-c", "cor")) ){
    	if(dataType(obj) == "binary")
    		test.type = test.type[1]
    	else
    		test.type = test.type[2]
    }
    
	# extract the parameters used to call the function
	params = list()
	for(p in names(formals())){
		if( p == "" | p == "..." | p == "obj")
			next
		
		params[[p]] = eval(parse(text=p))
	}
	
	# save all relevant stats to be included in the DDGraph object
	stats = list()
	
	if(verbose){	
		message("Using parameters:")
		message("  name = ", datasetName(obj))
		message("  alpha = ", alpha)
		message("  p.value.adjust.method = ", p.value.adjust.method)
		message("  test.type = ", test.type)
		message("  max.set.size = ", max.set.size)
		message("  mc.replicates =" , mc.replicates)
		message("  report.file = ", report.file)
		message("  star = ", star)
		message("  min.table.size = ", min.table.size, "\n")
		
	}
	
	# the whole dataframe with class labels
	dd = rawData(obj)
	# only the input variables (without class labels)
	vars = dd[,1:(ncol(dd)-1)]
	
	num.vars = ncol(vars)
	
	# calculate the maximal size of the conditioning set if applicable	
	if( dataType(obj) == "binary" & is.null(max.set.size) ){
		minority.sum = sum(dd$class)
		if(minority.sum > length(dd$class)/2)
			minority.sum = length(dd$class) - minority.sum
		max.set.size = trunc( log2(minority.sum) - 2 )
		if(max.set.size < 0)
			max.set.size = 0
		if(verbose)
			message("Estimating the maximal conditioning set size to ", max.set.size, "\n")
			
		if(max.set.size == 0)
			warning("Too few samples in the positive class, skipping the calculation of conditional independence tests.")
	}
	
	# default is all available
	if(is.null(max.set.size) & dataType(obj) == "continuous"){
		max.set.size = ncol(vars) - 1
	}
	
	# error checking and edge cases
	if(is.null(max.set.size)){
		stop("Please provide the maximal conditioning set size max.set.size parameter")
	}	

	# helper function for making square matrices with named dimensions	
	matrixWithNames = function(value, size, data.names){
		ret = matrix(value, ncol=size, nrow=size)
		colnames(ret) = data.names
		rownames(ret) = data.names
		ret
	}	
	
	# start an empty report file
	if( !is.null(report.file) ){
		options(width=150)
		sink(report.file)		
		sink()
	}	
	
	# store the conditioning sets and associated p values
	ciTests = list()
	for(i in 1:ncol(vars))
		ciTests[[i]] = list()
		
	###
	#  1) Find enriched variable set E: test for X ⊥ C | ∅ where X is any variable and C is the class labels.)
	#####
	setE.pvals = structure(rep(NA, ncol(vars)), names=names(vars))
	setE.log2FC = structure(rep(NA, ncol(vars)), names=names(vars))
	
	for(i in 1:ncol(vars)){
		res = ciTest(obj, i, "class", test.type=test.type, B=mc.replicates, min.table.size=min.table.size)
		setE.pvals[i] = res$pValue
		if( dataType(obj) == "binary" )
			setE.log2FC[i] = log2( mean(dd[dd$class==1,i]) / mean(dd[dd$class==0,i]) )
		else
			setE.log2FC[i] = NA
		
		ciTests[[i]] = append(ciTests[[i]], res) # store the results of the ciTest
	}
	
	# adjust p values 
	setE.pvals = p.adjust(setE.pvals, method=p.value.adjust.method)
	
	# find variable in set E	
	setE.inx = which(setE.pvals <= alpha)
	setE.names = names(setE.pvals)[setE.inx]
	
	stats[["setE.pvals"]] = setE.pvals
	stats[["setE.log2FC"]] = setE.log2FC
	stats[["setE.selected"]] = setE.inx	
	stats[["ciTests"]] = ciTests
	
	# make the final calls stuff
	final.calls = data.frame("name"=names(vars), "type"="", "explained.by"="", "explained.pval"=0, 
		"conditional.type"="", "conditional.on"="", "conditional.explained"="", "conditional.pval"="",
		"marginal.pval"=setE.pvals, "log2FC"=setE.log2FC, stringsAsFactors=F)
	rownames(final.calls) = NULL
	
	# all that are in setE.inx are direct, otherwise non-informative
	final.calls$type = c("non-informative", "direct")[((1:ncol(vars)) %in% setE.inx)+1]
	final.calls$explained.pval = NA	
		
	stats[["final.calls"]] = final.calls	
	
	# check if we have some enriched variables at all
	if( length(setE.inx) == 0 ){
		return(new("DDGraph", dataset=obj, params=params, stats=stats, direct=NULL, 
                   indirect=NULL, joint=NULL, edges=list()))
	}
	# check if only one variable is enriched, in that case return DDGraph with it alone
	if( length(setE.inx) == 1 & !star ){
		return(new("DDGraph", dataset=obj, params=params, stats=stats, direct=setE.inx, 
                   indirect=NULL, joint=NULL, edges=list()))
	}
	
	###
	#  2.1) Find the neighbourhood (Adj) of X,Y ∊ E: test X ⊥ Y | ∅ and X ⊥  C | Y (and Y is enriched)
	#####
	adjX.p.vals = matrixWithNames(NA, size=ncol(vars), data.names=names(vars))
	adjX.star.p.vals = matrixWithNames(NA, size=ncol(vars), data.names=names(vars))
	adjX.star.cond.p.vals = matrixWithNames(NA, size=ncol(vars), data.names=names(vars))
	adjX.tests = list()
	# list of conditional ciTests, in two level-structure cond.ciTest[[i]][[j]] for test I(i, "class" | j)
	cond.ciTests = vector(mode="list", length=ncol(vars)) 
	for(i in 1:ncol(vars))
		cond.ciTests[[i]] = vector(mode="list", length=ncol(vars)) 

	for(i in 1:ncol(vars)){
		for(j in 1:ncol(vars)){
			# only the upper triangle, since  adjX.p.vals[i,j]==adjX.p.vals[j,i]
			if( j > i ){
				# DDGraph algorithm
				if( i %in% setE.inx && j %in% setE.inx){ 
					res = ciTest(obj, i, j, test.type=test.type, B=mc.replicates, min.table.size=min.table.size)
					adjX.p.vals[i,j] = res$pValue
					adjX.tests[[paste(i, j)]] = res
					if( star ){
						adjX.star.p.vals[i,j] = res$pValue
					}
				
				} else if( star ){ # DDGraph* algorithm
					res = ciTest(obj, i, j, test.type=test.type, B=mc.replicates, min.table.size=min.table.size)
					adjX.star.p.vals[i,j] = res$pValue					
				}
			}
			if(star){				
				# do I(i, "class" | j)
				if(i != j & j %in% setE.inx){
					cond.ciTests[[i]][[j]] = ciTest(obj, i, "class", j, test.type=test.type, B=mc.replicates, min.table.size=min.table.size)
					adjX.star.cond.p.vals[i,j] = cond.ciTests[[i]][[j]]$pValue
				}
			}
		}
	}
	
	# apply p value adjustment for the X adjacency matrix
	adjX.p.vals = p.adjust(adjX.p.vals, method=p.value.adjust.method)
	if( is.vector(adjX.p.vals) ){
		# we are using R >= 2.13, need to convert back to matrix		
		adjX.p.vals = matrix(adjX.p.vals, ncol=ncol(vars))
		colnames(adjX.p.vals) = names(vars)
		rownames(adjX.p.vals) = names(vars)
	}
	
	adjX.star.p.vals = p.adjust(adjX.star.p.vals, method=p.value.adjust.method)
	if( is.vector(adjX.star.p.vals) ){
		# we are using R >= 2.13, need to convert back to matrix		
		adjX.star.p.vals = matrix(adjX.star.p.vals, ncol=ncol(vars))
		colnames(adjX.star.p.vals) = names(vars)
		rownames(adjX.star.p.vals) = names(vars)
	}
	
	adjX.star.cond.p.vals = p.adjust(adjX.star.cond.p.vals, method=p.value.adjust.method)
	if( is.vector(adjX.star.cond.p.vals) ){
		# we are using R >= 2.13, need to convert back to matrix		
		adjX.star.cond.p.vals = matrix(adjX.star.cond.p.vals, ncol=ncol(vars))
		colnames(adjX.star.cond.p.vals) = names(vars)
		rownames(adjX.star.cond.p.vals) = names(vars)
	}
	
	# the neighbourhood of X
	adjX.inx = vector("list", ncol(vars))
	# set of nodes S for X: for y in S: I(x, "class"|y) = false (i.e. x is dependent on class conditioned on y)
	condX.inx = vector("list", ncol(vars))
	if(star){
		for(i in 1:ncol(vars)){
			for(j in 1:ncol(vars)){ 
				# only if !(i ⊥   j)
				if( j > i & adjX.star.p.vals[i,j] <= alpha){
					# only add if the !(j ⊥  class | i)
					if( !is.na(adjX.star.cond.p.vals[j,i]) && adjX.star.cond.p.vals[j,i] <= alpha )
						adjX.inx[[i]] = union(adjX.inx[[i]], j)
					# only add if the !(i ⊥  class | j)
					if( !is.na(adjX.star.cond.p.vals[i,j]) && adjX.star.cond.p.vals[i,j] <= alpha )	
						adjX.inx[[j]] = union(adjX.inx[[j]], i)
				}
				
				# the conditioning neighbhood
				if( !is.na(adjX.star.cond.p.vals[i,j]) && adjX.star.cond.p.vals[i,j] <= alpha )
					condX.inx[[i]] = union( condX.inx[[i]], j )
			}
		}
	}	
	
	# all variables that are associated with class conditional on another variable (now including setE)
	setE.cond = which( ! sapply(condX.inx, is.null) )  #setdiff(which( ! sapply(condX.inx, is.null) ), setE.inx)
	onlySetE.cond = setdiff(setE.cond, setE.inx) #  variables that are only conditional, but not primarly enriched
	names(setE.cond) = names(vars)[setE.cond]
	
	# reset the ciTests for conditional E set
	#if(length(setE.cond) > 0)
	#	ciTests[setE.cond] = vector("list", length(setE.cond))
	
	###
	#  2.2) Test if association between X,Y ∊ E can be explained by C: test X ⊥ Y | C
	#####
	expC.p.vals = matrixWithNames(NA, size=ncol(vars), data.names=names(vars))
	expC.tests = list()

	for(i in setE.inx){
		for(j in setE.inx){
			if( j > i ){ # only the upper triangle, since  adjX.p.vals[i,j]==adjX.p.vals[j,i]
				res = ciTest(obj, i, j, "class", test.type=test.type, B=mc.replicates, min.table.size=min.table.size)
				expC.p.vals[i,j] = res$pValue
				expC.tests[[paste(i,j)]] = res
			}
		}
	}
	
	# apply p value adjustment for the X adjacency matrix
	expC.p.vals = p.adjust(expC.p.vals, method=p.value.adjust.method)
	if( is.vector(expC.p.vals) ){
		# we are using R >= 2.13, need to convert back to matrix
		expC.p.vals = matrix(expC.p.vals, ncol=ncol(vars))
		colnames(expC.p.vals) = names(vars)
		rownames(expC.p.vals) = names(vars)
	}
	
	# report what we calculated so far
	if( !is.null(report.file) ){
		sink(report.file, append=T)
		cat("setE.pvals\n")
		print(setE.pvals)
		cat("\nsetE.inx\n")
		print(setE.inx)
		cat("\nsetE.names\n")
		print(setE.names)
		cat("\nadjX.p.vals\n")
		print(adjX.p.vals)
		cat("\nexpC.p.vals\n")
		print(expC.p.vals)
		if(star){
			cat("\nadjX.star.p.vals\n")
			print(adjX.star.p.vals)
			cat("\nadjX.star.cond.p.vals\n")
			print(adjX.star.cond.p.vals)
			cat("\nadjX.inx\n")
			print(adjX.inx)
			cat("\ncondX.inx\n")
			print(condX.inx)
			cat("\nsetE.cond\n")
			print(setE.cond)
			cat("\nonlySetE.cond\n")
			print(onlySetE.cond)
			cat("\ncond.ciTests\n")
			print(cond.ciTests)
		}
		sink()
	}
	
	###
	#   3) Start the main loop:
	#####
	
	# initialise the neighbourhood of C to the whole set E, contains indexes of variables adjecent to class labels (C)
	# this contains cond indep tests
	#adjC = union(setE.inx, setE.cond)
	adjC = list()
	for(i in setE.inx){
		adjC[[length(adjC)+1]] = ciTests[[i]][[1]]
	}
	for(i in setE.cond){
		for(j in condX.inx[[i]])
			adjC[[length(adjC)+1]] = cond.ciTests[[i]][[j]]
	}
	
	if( max.set.size == 0){		
		return(new("DDGraph", dataset=obj, params=params, stats=stats, direct=setE.inx, 
                   indirect=NULL, joint=NULL, edges=list()))
	}
	
	# update this to the estimated value if any
	params[["max.set.size"]] = max.set.size 
	
	#####
	#    3.1) Set n = 1 (size of conditioning set), 	
	#    3.2) For each X ∊ Adj(C) calculate X ⊥ C | S where |S|=n and S ⊆ Adj(C)\X or S ⊆ Adj(X) (DDGraph*)
	#       find the largest p.value for all values of S and test if p.value > alpha, if so, 
	#       record Dsep(X,C,S) and remove X from Adj(C)
	#    3.3) increment n=n+1 until n = log_2 (#positives) - 2 or n=|adj(C)|
	#############
	
	# store the indexes to those sets S which d-separate best this each X from C
	Dsep = structure(vector(mode="list", length=ncol(vars)), names=names(vars)) 
	adjC.pvals.at.n = list() # the adjC p values at different values of conditioning set size (n) - BEFORE multiple testing correction
	adjC.at.n = list() # the adjC neighbouhood at different conditioning set sizes (n)
	adjC.pvals.adj.at.n = list() # the adjusted adjC p values at different n values
	ciTests.adjC.removed = list() # list of removed ciTests.adjC elements 
	adjC.initial = adjC # the initial value of adjC we started with
	
	# initiate all the ciTest for adjC, id -> tests
	ciTests.adjC = list()
	for(x in adjC)
		ciTests.adjC[[ CITestResultID(x) ]] = list(x)
		
	# list of joint conditional independence tests
	jointTests = list()
	
	for( n in 1:max.set.size ){
		# only one in adjacency, stop
		if(length(adjC) <= 1)
			break		
		
		# break if we reduced the adjacency below the size of conditioning set	
		if(n >= adjC.condSetSize(adjC) )
			break
			
		if( !is.null(report.file) ){
			sink(report.file, append=T)
			cat("\nConditioning on sets of", n, "variable(s) from", paste(adjC.allVarNames(adjC), collapse=", "), "\n")
			sink()
		}
		
		# initialize S.comb.adjC
		if(star){
			# make all combination sizes up to n, this is only for star, otherwise in unneccessary
			S.comb.adjC = list()
			for(i in 1:n){
				l = as.list(as.data.frame(t(combinations(length(adjC), i))))
				names(l) = NULL
				S.comb.adjC = append(S.comb.adjC, l)
			}
		} else { 
			l = as.list(as.data.frame(t(combinations(length(adjC), n))))
			names(l) = NULL
			S.comb.adjC = l
		}
		
		# x is an CITestResult object
		for(x in adjC){
			x.inx = x@targetInx
			x.cond = x@condSetInx
			# generate combinations of nodes adjacent to C (of size n)		
			for(i in 1:length(S.comb.adjC)){
				S = adjC [ S.comb.adjC[[i]] ]
				S.inx = adjC.allVarInx(S)
				
				# skip subsets with X in them..
				if(x.inx %in% S.inx )  
					next

				S.inx = union(S.inx, x.cond)

				# check if we already did this conditional independence test				
				if(star){
					done.already = sapply(ciTests.adjC[[CITestResultID(x)]], function(y) setequal(y@condSetInx, S.inx))
					if(any(done.already))
						next
				}
					
				# just process those conditional set sizes that are equal to n
				if(length(S.inx) != n)
					next 
				
				# check if we already done this test before for another conditioning set
				if(star)
					done.already = sapply(ciTests[[x.inx]], function(y) setequal(y@condSetInx, S.inx))
				else
					done.already = FALSE
					
				if(any(done.already)){
					stopifnot(sum(done.already) == 1)
					res = ciTests[[x.inx]][[which(done.already)]]
				} else {
					# doing the test for the first time, record
					res = ciTest(obj, x.inx, "class", S.inx, test.type=test.type, B=mc.replicates, min.table.size=min.table.size)
					ciTests[[x.inx]] = append(ciTests[[x.inx]], res)
				}
				
				ciTests.adjC[[ CITestResultID(x) ]] = append(ciTests.adjC[[ CITestResultID(x) ]], res)
								
			}		
		}
		
		# extract the p values from the CI tests
		pvals.list = extractCITestResultProperty(ciTests.adjC, "pValue")
		# pvals.list = pvals.list[adjC]
		
		# for each variable in adjC find the maximal p value
		adjC.pvals = unlist(lapply( pvals.list, max ))
		adjC.pvals.inx = unlist(lapply( pvals.list, which.max ))
		
		# note that these vectors are not of size ncols(vars) but of length(adjC)
		#stopifnot(length(adjC.pvals) == length(adjC))
		#stopifnot(length(adjC.pvals.inx) == length(adjC))
		
		# record the adjC values before multiple test correction
		adjC.at.n[[n]] = adjC
		adjC.pvals.at.n[[n]] = adjC.pvals
		
		if( !is.null(report.file) ){
			sink(report.file, append=T)
			cat("\nadjC at n =", n, "\n")
			print(adjC)
			cat("\nciTests.adjC at n =", n, "\n")
			print(ciTests.adjC)
			cat("\nadjC.pvals at n =", n, "\n")
			print(adjC.pvals)
			cat("\nadjC.pvals.inx at n =", n, "\n")
			print(adjC.pvals.inx)
			#for(ii in adjC){
			#	print(ciTests[[ ii ]][[adjC.pvals.inx[ii] ]])
			#}
			sink()
		}
		
		# do multiple test correction if any
		adjC.pvals = p.adjust(adjC.pvals, method=p.value.adjust.method)
		
		if( !is.null(report.file) ){
			sink(report.file, append=T)
			cat("\nadjC.pvals (adjusted) at n =", n, "\n")
			print(adjC.pvals)
			sink()
		}
		
		adjC.pvals.adj.at.n[[n]] = adjC.pvals
		
		# delete edges that went over the p value treshold
		adjC.delete = which(adjC.pvals > alpha)
		to.delete = rep(FALSE, length(adjC))
		delete.inx = c()
		
		stopifnot(length(adjC.pvals) == length(adjC))
		
		# delete stuff from adjC
		if( length(adjC.delete) > 0 ){			
			to.delete[adjC.delete] = TRUE
			
			# collect deleted citests			
			delete.citests = list()
			delete.citests.all = list()
			for(i in adjC.delete){
				delete.citests[[ length(delete.citests)+1 ]] = ciTests.adjC[[i]][[ adjC.pvals.inx[i] ]]
				
				# find those tests that when substituted have P-value larger than threshold				
				above.alpha = sapply(ciTests.adjC[[i]], function(indtest){ 
					## see if it would be deleted if substituted 
					adjC.pvals.del = adjC.pvals
				    adjC.pvals.del[i] = indtest@pValue
				    
					adjC.pvals.del = p.adjust(adjC.pvals.del, method=p.value.adjust.method)

					#indtest@pValue > alpha 
					adjC.pvals.del[i] > alpha
					}
				)
				
				# change the P-value to the last corrected P-value if its above alpha
				# since it is going to be removed
				for(aa in which(above.alpha)){
					indtest = ciTests.adjC[[i]][[aa]]
					
					# corrected P-value
					adjC.pvals.del = adjC.pvals
				    adjC.pvals.del[i] = indtest@pValue				    
					adjC.pvals.del = p.adjust(adjC.pvals.del, method=p.value.adjust.method)
					
					
					indtest@pValue = adjC.pvals.del[i]
					ciTests.adjC[[i]][[aa]] = indtest
					
				}
				
				delete.citests.all = c(delete.citests.all, ciTests.adjC[[i]][ above.alpha ])
			}
			
			# detected joint tests
			# need to use ALL of them, not only the max ones
			for(i in 1:length(delete.citests.all)){
				for(j in 1:length(delete.citests.all)){
					if(i > j && CITestResultID(delete.citests.all[[i]]) != CITestResultID(delete.citests.all[[j]])){
						# they are joint if the set of targets and condsetinx is the same
						set1 = c(delete.citests.all[[i]]@targetInx, delete.citests.all[[i]]@condSetInx)
						set2 = c(delete.citests.all[[j]]@targetInx, delete.citests.all[[j]]@condSetInx)
					
						#stopifnot(length(delete.citests.all[[i]]@condSetInx) == length(delete.citests.all[[j]]@condSetInx))
					
						if(setequal(set1, set2)){
							jointTests[[length(jointTests)+1]] = list(delete.citests.all[[i]], delete.citests.all[[j]])
						}
					}
				}
			}
			
			
			# collect all unconditional variables to be deleted and register Dseps
			delete.inx = c()
			for(i in 1:length(adjC)){
				x = adjC[[i]]
				if(i %in% adjC.delete && length(x@condSetInx) == 0){
					delete.inx = c(delete.inx, x@targetInx)
					
					# get the citest responsible for deletion
					citest = ciTests.adjC[[i]][[ adjC.pvals.inx[i] ]]
					# identify this citest in ciTests
					real.x.inx = citest@targetInx
					Dsep[[ real.x.inx ]] = citest
				}
			}
			
			stopifnot( setequal(delete.inx, unique(delete.inx)) )
			
			# remove all conditional that are conditionong in those variables removed
			for(i in 1:length(adjC)){
				x = adjC[[i]]
				if(length(x@condSetInx)>0 && all(x@condSetInx %in% delete.inx))
					to.delete[i] = TRUE
			}
			
			for(i in which(to.delete)){
				citests = ciTests.adjC[[i]]
				key = CITestResultID(adjC[[i]]) # ID
				ciTests.adjC.removed[[ key ]] = citests
			}
					
			adjC = adjC[!to.delete]
			ciTests.adjC = ciTests.adjC[!to.delete]
		}
		
		if( !is.null(report.file) ){
			sink(report.file, append=T)
			cat("\nadjC.delete\n")
			print(adjC.delete)
			cat("\ndelete.inx\n")
			print(delete.inx)
			cat("\nto.delete\n")
			print(to.delete)
			cat("\njointTests\n")
			print(jointTests)
			cat("\nciTests.adjC.removed\n")
			print(ciTests.adjC.removed)
			cat("\nVariables left after iteration", n, ":", paste(adjC.allVarNames(adjC), collapse=", "), "\n")
			cat("\nadjC\n")
			print(adjC)
			cat("======================================\n")
			sink()
		}
	}
	
	# end of main loop
################################################################################################################
################################################################################################################
	
	## remove those from adjC that are both direct and conditional
	adjC.uncond = c()
	for(x in adjC){
		if(length(x@condSetInx) == 0)
			adjC.uncond = c(adjC.uncond, x@targetInx)
	}
	# remove conditional on those adjC.uncond
	adjC.filtered = list()
	for(x in adjC){
		if(!(x@targetInx %in% adjC.uncond & length(x@condSetInx) > 0)){
			adjC.filtered = c(adjC.filtered, x)
		}
	}
	adjC.before.filtering = adjC
	adjC = adjC.filtered
	
	## similarly filter out joint tests
	jointTests.filtered=list()
	for(x in jointTests){
		if( !(x[[1]]@targetInx %in% adjC.uncond) && !(x[[2]]@targetInx %in% adjC.uncond) ){
			jointTests.filtered[[length(jointTests.filtered)+1]] = x
		}
	}
	jointTests.before.filtering = jointTests
	jointTests = jointTests.filtered
	
	## find joint variables - only among the indirect
	direct.nodes = c()
	indirect.nodes = c()
	joint.nodes = c()
	conditional.nodes = c()
	conditional.joint.nodes = c()
	# index all the tests for each of the joint variables
	jointNodeTests = list()
	condJointNodeTests = list()
	
	# find the direct and conditional
	for(x in adjC){
		node.inx = x@targetInx
		names(node.inx) = x@targetName
		if(length(x@condSetInx) == 0){
			direct.nodes = c(direct.nodes, x)
		} else {
			conditional.nodes = c(conditional.nodes, x)			
		}
	}
	
	# multiple testing correction might remove a node at later stage
	# make sure we catch this and kick out the node
	to.remove.from.indirect = c()
	if(p.value.adjust.method != "none" & length(ciTests.adjC.removed)>0){		
		to.remove.from.adjC.removed = c()
		for(i in 1:length(ciTests.adjC.removed)){
			citests = ciTests.adjC.removed[[i]]
			for(j in 1:length(citests)){
				if(length(citests[[j]]@condSetInx) == 0 && citests[[j]]@pValue > alpha){
					to.remove.from.indirect = c(to.remove.from.indirect, citests[[j]]@targetInx)
					to.remove.from.adjC.removed = c(to.remove.from.adjC.removed, i)
				}
			}
		}	
		
		if(length(to.remove.from.adjC.removed)>0){
			ciTests.adjC.removed = ciTests.adjC.removed[-to.remove.from.adjC.removed]
			# loop throught joint tests and do the same
			to.remove.joint.tests = c()
			for(j in 1:length(jointTests)){
				if((jointTests[[j]][[1]]@targetInx %in% to.remove.from.indirect) | 
				   (jointTests[[j]][[2]]@targetInx %in% to.remove.from.indirect)){
					to.remove.joint.tests = c(to.remove.joint.tests, j)
				}
			}
			if(length(to.remove.joint.tests) > 0)
				jointTests = jointTests[-to.remove.joint.tests]
				
			# also remove all conditioning not to confuse the algorithm later on.. 
			for(j in 1:length(ciTests.adjC.removed)){
				to.remove.cond = which(sapply(ciTests.adjC.removed[[j]], function(x) to.remove.from.indirect %in% x@condSetInx))
				ciTests.adjC.removed[[j]] = ciTests.adjC.removed[[j]][-to.remove.cond]				
			}
		}
	}
	
	if( !is.null(report.file) ){
		sink(report.file, append=T)
		cat("\nciTests\n\n")
		print(ciTests)
		cat("\nadjC.uncond\n\n")
		print(adjC.uncond)
		cat("\nadjC.filtered\n\n")
		print(adjC.filtered)
		cat("\nadjC.before.filtering\n\n")
		print(adjC.before.filtering)
		cat("\njointTests.before.filtering\n")
		print(jointTests.before.filtering)
		cat("\njointTests\n\n")
		print(jointTests)
		sink()
	}
	
	# removed - either indirect or joint
	removed.ids = names(ciTests.adjC.removed)
	
	direct.inx = extract.targetInx(direct.nodes)
	
	##########################################################
	# find the indirect and joint
	########################################################
	all.joint.ids = adjC.toIDs(unlist(jointTests))
	for(x in removed.ids){
		stopifnot(!is.null(ciTests.adjC.removed[[x]]))
		
		# see if all conditional independencies are joint or not
		citests = ciTests.adjC.removed[[x]]
		above.alpha = sapply(citests, function(xx) xx@pValue > alpha)
		ids = adjC.toIDs(citests[above.alpha])
		
		if(length(ids)>0 && length(all.joint.ids)>0 && all(ids %in% all.joint.ids)){
			# decide if it is conditionally joint or not
			if(length(citests[[1]]@condSetInx)>0){
				# if the conditoning node is not direct, don't register as cond joint
				if(all(citests[[1]]@condSetInx %in% direct.inx)){
					conditional.joint.nodes = c(conditional.joint.nodes, citests[[1]])
					condJointNodeTests[[length(condJointNodeTests)+1]] = citests[above.alpha]
				}
			} else {
				joint.nodes = c(joint.nodes, citests[[1]])
				jointNodeTests[[length(jointNodeTests)+1]] = citests[above.alpha]
			}
		} else {
			if(length(ciTests.adjC.removed[[x]][[1]]@condSetInx)==0)
				indirect.nodes = c(indirect.nodes, ciTests.adjC.removed[[x]][[1]])
		}
	}
	
	if(verbose){
		message("Direct: ", paste(sapply(direct.nodes, CITestResultVar), collapse=" "))
		message("Joint: ", paste(sapply(joint.nodes, CITestResultVar), collapse=" "))
		message("Indirect: ", paste(sapply(indirect.nodes, CITestResultVar), collapse=" "))
		message("Conditional: ", paste(sapply(conditional.nodes, CITestResultVar), collapse=" "))
		message("Conditional Joint: ", paste(sapply(conditional.joint.nodes, CITestResultVar), collapse=" "))
	}
	
	# indicies of direct variables
		
	joint.inx = extract.targetInx(joint.nodes)
	indirect.inx = extract.targetInx(indirect.nodes)
	conditional.inx = extract.targetInx(conditional.nodes)
	conditional.joint.inx = extract.targetInx(conditional.joint.nodes)
	
	#######################
	# Construct a DDGraph, start by taking new citest to show for indirect variables
	#
	###########################
	
	#######################
	#  Find indirect variables Dseps to show
	##############################	
	for(x in indirect.nodes){
		# try to take those that are direct
		citests = ciTests.adjC.removed[[ CITestResultID(x) ]]
		stopifnot(!is.null( citests ))
		
		# make the direct/joint for Dsep
		if(length(citests[[1]]@condSetInx)==0){
			x.inx = citests[[1]]@targetInx		
			direct.and.above.cutoff = rep(FALSE, length(citests))
			joint.and.above.cutoff = rep(FALSE, length(citests))
			for(i in 1:length(citests)){
				if(all(citests[[i]]@condSetInx %in% direct.inx) & citests[[i]]$pValue > alpha)
					direct.and.above.cutoff[i] = TRUE
				if(all(citests[[i]]@condSetInx %in% joint.inx) & citests[[i]]$pValue > alpha)
					joint.and.above.cutoff[i] = TRUE
			}
			
			ids = adjC.toIDs(citests)
			ids.joint = ids %in% all.joint.ids
			
			# make sure we don't pick up any of the joint relationships
			if(any(ids.joint))
				joint.and.above.cutoff[ids.joint] = FALSE
		
			if(any(direct.and.above.cutoff)){
				p.inx = which(direct.and.above.cutoff)
				pvals = sapply(citests[p.inx], function(x) x@pValue)
				max.inx = p.inx[which.max(pvals)]
				
				Dsep[[x.inx]] = citests[[max.inx]]
			} else if(any(joint.and.above.cutoff)){
				p.inx = which(joint.and.above.cutoff)
				pvals = sapply(citests[p.inx], function(x) x@pValue)
				max.inx = p.inx[which.max(pvals)]
				
				Dsep[[x.inx]] = citests[[max.inx]]

			} else {
				# select any non-joint for Dsep
				alpha.rank = rank(sapply(citests, function(xx) xx@pValue))
				
				alpha.rank[ids.joint] = 0
				max.inx = which.max(alpha.rank)
				#if(citests[[max.inx]]@pValue <= alpha){
				#	browser()
				#}
				stopifnot(citests[[max.inx]]@pValue > alpha)
				
				Dsep[[x.inx]] = citests[[max.inx]]
			}
		}		
	}
	
	## make a report on Dseps
	if(verbose){
		cat("Best explaining conditional independencies:\n")
		for(i in 1:length(Dsep)){
			if(!is.null(Dsep[[i]])){
				dsep = Dsep[[i]]
				cat(dsep@targetName, "<-", paste(dsep@condSetName, collapse=","), "with pval =", dsep@pValue, "\n")
			}
		}
	}
	
	if( !is.null(report.file) ){
		sink(report.file, append=T)
		cat("\nadjC\n")
		print(adjC)
		cat("\nremoved.ids\n")
		print(removed.ids)
		cat("\ndirect.nodes\n")
		print(direct.nodes)
		cat("\nindirect.nodes\n")
		print(indirect.nodes)
		cat("\njoint.nodes\n")
		print(joint.nodes)
		cat("\nconditional.nodes\n")
		print(conditional.nodes)
		cat("\ndirect.inx\n")
		print(direct.inx)
		cat("\nindirect.inx\n")
		print(indirect.inx)
		cat("\njoint.inx\n")
		print(joint.inx)
		cat("\nconditional.inx\n")
		print(conditional.inx)
		cat("\nDsep\n")
		print(Dsep)
		cat("\nall.joint.ids\n")
		print(all.joint.ids)
		sink()
	}
	
	##### Make DDGraph edges and such
	
	##########
	# Construct the edge set for the DDGraph
	#
	# 6) Construct an e-graph as follows:
	#    6.1) classify the variables from E as follows:
	#      - X is direct if X ∊ Adj(C)
	#      - X is indirect if Dsep(X,C,S) and not Joint(X,S)
	#      - X is joint otherwise
	#    6.2) draw a directed edge S -> X if Dsep(X,C,S) and not Joint(X,S)
	#    6.3) draw a bidirectonal edge S <-> X if Joint(X,S)
	#    6.4) for every X,Y ∊ adj(C) draw:
	#         - no edge X  Y if X ⊥ Y | ∅ (i.e. X not in Adj(Y))
	#         - dashed edge X - - Y if X ⊥ Y | C
	#         - undirected full edge X - Y otherwise
	####################


	# construct a list of edges
	edges = list()
	# this will contain all edges, even those that won't be drawn (i.e. when directed takes precedence over joint)
	allEdges = list() 
	edgesDrawn = matrixWithNames("none", size=ncol(vars), data.names=names(vars))
		
	# 6.2) add all directional edges in Dsep
	for(x in indirect.nodes){
		dsep = Dsep[[ x@targetInx ]]
		if(length(dsep$condSetInx)>0){
			# check if this Dsep is part of a joint connection, if so, don't draw it as directed line
			if( !(CITestResultID(dsep) %in% all.joint.ids) ){
				for(i in 1:length(dsep$condSetInx)){
					edge = new("DDGraphEdge", fromInx=dsep$condSetInx[i], fromName=dsep$condSetName[i],
					  toInx=dsep$targetInx, toName=dsep$targetName, ciTests=list(dsep), type="directed")
					
					edgesDrawn[edge@fromInx, edge@toInx] = "directed"
					edges = append(edges, edge)					
					allEdges = append(allEdges, edge)
				}
			}
		}
	}	
	
	# 6.3) add all the joint edges	
	if(length(jointTests) > 0){
		for(i in 1:length(jointTests)){
			st = jointTests[[i]][[1]]
			ot = jointTests[[i]][[2]]
			if(st$targetInx %in% c(joint.inx, conditional.joint.inx) | 
				ot$targetInx %in% c(joint.inx, conditional.joint.inx)){			
				edge = new("DDGraphEdge", fromInx=st$targetInx, fromName=st$targetName,
						   toInx=ot$targetInx, toName=ot$targetName, ciTests=list("joint"=st, "original"=ot), 
					       type="bidirectional")

				# avoid duplicates
				if(edgesDrawn[edge@fromInx, edge@toInx] != "bidirectional" 
					& edgesDrawn[edge@toInx, edge@fromInx] != "bidirectional"
					& edgesDrawn[edge@fromInx, edge@toInx] == "none" & edgesDrawn[edge@toInx, edge@fromInx] == "none"){
					edgesDrawn[edge@fromInx, edge@toInx] = "bidirectional"
					edgesDrawn[edge@toInx, edge@fromInx] = "bidirectional"
					edges = append(edges, edge)
				}
			
				allEdges = append(allEdges, edge)
			}
		}
	}	
	
	
	# 6.4) for every X,Y ∊ (direct union joint) draw: no edge, dashed edge, undirected edge
	for(x in union(direct.inx, joint.inx)){
		for(y in union(direct.inx, joint.inx)){
			# only for upper triangle since these are all symmetric relationships
			if( y > x & edgesDrawn[x,y] == "none" & edgesDrawn[y,x] == "none"){
				if( adjX.p.vals[x,y] >= alpha ){
					# should remain no edge
					edgesDrawn[x,y] = "none" 
				} else if( expC.p.vals[x,y] >= alpha ){
					# dashed					
					edgesDrawn[x,y] = "dashed"
					edgesDrawn[y,x] = "dashed"
					key = paste(x,y)
					edge = new("DDGraphEdge", fromInx=x, fromName=names(vars)[x],
					  toInx=y, toName=names(vars)[y], ciTests=list(adjX.tests[[key]], expC.tests[[key]]), 
					  type="dashed")
					edges = append(edges, edge)
					allEdges = append(allEdges, edge)
				} else {
					# undirected
					edgesDrawn[x,y] = "undirected"
					edgesDrawn[y,x] = "undirected"
					key = paste(x,y)
					edge = new("DDGraphEdge", fromInx=x, fromName=names(vars)[x],
					  toInx=y, toName=names(vars)[y], ciTests=list(adjX.tests[[key]], expC.tests[[key]]), 
					  type="undirected")
					edges = append(edges, edge)
					allEdges = append(allEdges, edge)
				}
					
			}
		}
	}
	
	# 6.4a) add all conditional edges, connect the appropriate direct/joint with conditional
	for(x in conditional.nodes){
		x.inx = x@targetInx
		ys = x@condSetInx
		if(length(ys) > 0){
			for(y in ys){
				edgesDrawn[x.inx,y] = "undirected"
				edgesDrawn[y,x.inx] = "undirected"
				edge = new("DDGraphEdge", fromInx=x@targetInx, fromName=x@targetName,
					  toInx=y, toName=names(vars)[y], ciTests=list(x), 
					  type="dotted")
				edges = append(edges, edge)
				allEdges = append(allEdges, edge)
			}
		
		}
	}
	for(x in conditional.joint.nodes){
		x.inx = x@targetInx
		ys = x@condSetInx
		if(length(ys) > 0){
			for(y in ys){
				if(!(y %in% direct.inx))
					next
				edgesDrawn[x.inx,y] = "undirected"
				edgesDrawn[y,x.inx] = "undirected"
				edge = new("DDGraphEdge", fromInx=x@targetInx, fromName=x@targetName,
					  toInx=y, toName=names(vars)[y], ciTests=list(x), 
					  type="dotted")
				edges = append(edges, edge)
				allEdges = append(allEdges, edge)
			}
		
		}
	}
	
	
	## provide summaries about edges in stats
	for(i in 1:nrow(final.calls)){
		if(! (i %in% setE.inx) | i %in% to.remove.from.indirect){
			final.calls$type[i] = "no dependence"
			final.calls$explained.pval[i] = NA
		} else if(i %in% direct.inx){
			final.calls$type[i] = "direct"
			final.calls$explained.pval[i] = NA
		} else if(i %in% indirect.inx){
			final.calls$type[i] = "indirect"
			inx = which( sapply(allEdges, function(x) x@toInx==i & x@type=="directed") )[1]
			# the original conditional independence test
			if( is.null(allEdges[[inx]]) ){
				cat("Internal ERROR with finding the edge for indirect variable\n")
				print(i)
				print(allEdges)
				stop("Internal ERROR with finding the edge for indirect variable")
			} 
			res = allEdges[[inx]]@ciTests[[1]] 
			final.calls$explained.by[i] = paste(res$condSetName, collapse=", ")
			final.calls$explained.pval[i] = res$pValue
			if( all(res$condSetInx %in% direct.inx) )
				final.calls$type[i] = "strong indirect"
			else
				final.calls$type[i] = "weak indirect"
		} else {
			final.calls$type[i] = "joint"
			res = jointNodeTests[[ which(joint.inx == i) ]][[1]]
			final.calls$explained.by[i] = paste(res$condSetName, collapse=", ")
			final.calls$explained.pval[i] = res$pValue
		}
		
		# work out the conditional type
		if(i %in% conditional.inx){
			cond.inx = which(i == sapply(conditional.nodes, function(x) x@targetInx))
			final.calls$conditional.type[i] = paste(rep("direct", length(cond.inx)), collapse=" | ")
			for(dsep in conditional.nodes[ cond.inx ] ){
				ys = dsep@condSetInx
				if(final.calls$conditional.on[i] == "")
					final.calls$conditional.on[i] = CITestResultVar(dsep)
				else
					final.calls$conditional.on[i] = paste(final.calls$conditional.on[i], 
														  CITestResultVar(dsep), sep=" | ")
			}
		} 
		if(i %in% conditional.joint.inx){
			num.cond.joint = length(which(conditional.joint.inx == i))
			if(final.calls$conditional.type[i] != "")
				final.calls$conditional.type[i] = paste( final.calls$conditional.type[i],
					paste(rep("joint", num.cond.joint), collapse=" | "), sep=" | " )
			else
				final.calls$conditional.type[i] = paste(rep("joint", num.cond.joint), collapse=" | ")
			
			for(joint.match in which(conditional.joint.inx == i)){
				dsep = conditional.joint.nodes[[ joint.match ]]				
				#final.calls$explained.by[i] = paste(dsep$condSetName, collapse=", ")
				#final.calls$explained.pval[i] = dsep$pValue
				if(final.calls$conditional.on[i] == "")
					final.calls$conditional.on[i] = CITestResultVar(dsep)
				else
					final.calls$conditional.on[i] = paste(final.calls$conditional.on[i], 
															  CITestResultVar(dsep), sep=" | ")															  
															  
				# fill out explain and pval
				citests = condJointNodeTests[[ joint.match ]]
				if( final.calls$conditional.explained[i] == "")
					final.calls$conditional.explained[i] = paste(
						sapply(citests, function(x) paste(setdiff(x@condSetName, dsep@condSetName), collapse=",")),
						collapse=" | ")
				else
					final.calls$conditional.explained[i] = paste(final.calls$conditional.explained[i], paste(
						sapply(citests, function(x) paste(setdiff(x@condSetName, dsep@condSetName), collapse=",")),
						collapse=" | "), sep=" | ")
						
				if( final.calls$conditional.pval[i] == "")
					final.calls$conditional.pval[i] = paste(
						sapply(citests, function(x) x@pValue),collapse=" | ")	
				else
					final.calls$conditional.pval[i] = paste(final.calls$conditional.pval[i], 
						paste(sapply(citests, function(x) x@pValue),collapse=" | "), sep=" | ")
			}
		}		
		
	}
	rownames(final.calls) = NULL
	
	stats[["final.calls"]] = final.calls
	stats[["ciTests"]] = ciTests
	
	# make sure everything is properly names
	if(length(direct.inx)>0)
		names(direct.inx) = names(vars)[direct.inx]
	if(length(indirect.inx)>0)
		names(indirect.inx) = names(vars)[indirect.inx]
	if(length(joint.inx)>0)
		names(joint.inx) = names(vars)[joint.inx]
	if(length(conditional.inx)>0)
		names(conditional.inx) = names(vars)[conditional.inx]
	if(length(conditional.joint.inx)>0)
		names(conditional.joint.inx) = names(vars)[conditional.joint.inx]
	
	if( !is.null(report.file) ){
		sink(report.file, append=T)
		cat("\nedges\n")
		print(edges)
		cat("\nstats\n")
		print(stats)
		sink()
	}
	
	new("DDGraph", dataset=obj, edges=edges, params=params, stats=stats, 
	    direct=direct.inx, indirect=indirect.inx, joint=joint.inx,
	    conditional=conditional.inx, conditionalJoint=conditional.joint.inx)
	
}

#' This is a helper function for \code{DDDataSet::ncpc()}. From a list of ciTestResult object extract 
#' a list containing only one property. 
#'
#' @title Extract CITestResult properties
#'
#' @param ciTestList a two-level list of ciTestResult objects
#' @param prop.name the name of the property to extract (one of the slot names)
#'
#' @return a vector with the extracted property
extractCITestResultProperty = function(ciTestList, prop.name){	
	res = lapply( ciTestList, function(obj){ lapply(obj, function(x) x[[prop.name]])})
	
	# unlist only if it's not the multiple fields
	if( prop.name != "condSetName" & prop.name != "condSetInx")
		res = lapply( res, unlist)
	
	return(res)
}

#' This function is only for DDGraph with multiple testing correction enabled. The overall procedure is similar to
#' that described in (Li&Wang 2009). This is a helper function for \code{DDDataSet:ncpc()}. The single P-value
#' of D-separation is substituted in the list of P-values, P-values adjusted and the resulting P-value after correction
#' in the context of other P-values reported. 
#'
#' @title Multiple testing correction procedure for ncpc()
#'
#' @param dsep the conditional independence test result (of type \code{CITestResult})
#' @param x the index of the variables
#' @param adjC.pvals.at.n the p values associated with the variables at size n of conditioning set (list [[n]] -> [pvals])
#' @param p.value.adjust.method the p value adjustment method (same as in p.adjust())
#'
#' @return the p value after multiple test correction (if any)
#' @references J. Li and Z. J Wang, "Controlling the false discovery rate of the association/causality structure 
#'  learned with the PC algorithm" The Journal of Machine Learning Research 10 (2009): 475-514.
pValueAfterMultipleTesting = function(dsep, x, adjC.pvals.at.n, p.value.adjust.method){
	# nothing to do if there is no multiple test correction
	if(p.value.adjust.method == "none")
		return(dsep$pValue)
	
	# size of conditioning set			
	n.val = length(dsep$condSetInx) 

	# substitute the p value
	adjC.pvals = adjC.pvals.at.n[[n.val]]
	adjC.pvals[x] = dsep$pValue

	adjC.pvals = p.adjust(adjC.pvals, method=p.value.adjust.method)
	
	return(adjC.pvals[x])
}

#' Returns the total size of conditioning set for adjC (i.e. all variables present in adjC)
#'
#' @param adjC the adjC list of conditional independence tests for variables "adjacent" to target variable C
#' @return sum of all conditioning set sizes plus size of adjC, i.e. all variables present in adjC
adjC.condSetSize = function(adjC){
	length(adjC.allVarNames(adjC))
}

#' Get all the variable names in adjC, both target and condSet
#'
#' @param adjC the adjC list of conditional independence tests for variables "adjacent" to target variable C
#' @return character vector (unique names)
adjC.allVarNames = function(adjC){
	unique(unlist(sapply(adjC, function(x) c(x@targetName, x@condSetName))))
}

#' Get all the variable indicies in adjC, both target and condSet
#'
#' @param adjC the adjC list of conditional independence tests for variables "adjacent" to target variable C
#' @return numeric vector (unique values)
adjC.allVarInx = function(adjC){
	unique(unlist(sapply(adjC, function(x) c(x@targetInx, x@condSetInx))))
}

#' Get all the targetInx values in adjC
#'
#' @param adjC the adjC list of conditional independence tests for variables "adjacent" to target variable C
#' @return numeric vector (unique values)
adjC.targetInx = function(adjC){
	unique(unlist(sapply(adjC, function(x) x@targetInx)))
}

#' Make a list of conditional independence tests and converts them to IDs
#' @param adjC a list of conditional independence tests
adjC.toIDs = function(adjC){
	sapply(adjC, CITestResultID)
}

#' Extract all values of targetInx from a list of CITestResult
#'
#' @param adjC a list of CITestResult
extract.targetInx = function(adjC){
	ret = sapply(adjC, function(x) structure(x@targetInx, names=x@targetName))
	if(length(ret) == 0)
		ret = c()
		
	return(ret)
}


