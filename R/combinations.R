#
# Testing variable value combinations for binary data
#


#' Calculate which combinations of values of variables are significantly different in the 
#' two classes (only for binary data). This function takes an DDDataSet and a number of variables
#' and finds those combinations of values of those variables that have significantly different
#' frequencies in the two class labels. 
#'
#' @title Significant combinations of variables
#'
#' @param obj DDDataSet object
#' @param selected.vars indexes or names of variables selected for the test
#' @param cutoff the p value cutoff for reporting (default: 0.05)
#' @param p.adjust.method the multiple adjustment method (default: none)
#' @param verbose if to print progress output and additional information
#'
#' @return \code{data.frame} with ordered combinatorial patterns of selected variables
#' @export
#' @examples
#' data(mesoBin)
#' # find significant differences at 0.2 FDR
#' combinationsTest(mesoBin$Meso, c("Twi 2-4h", "Tin 6-8h", "Mef2 6-8h"), 0.2, "fdr")
combinationsTest = function(obj, selected.vars, cutoff=0.05, p.adjust.method="none", verbose=TRUE){
	if(dataType(obj) != "binary")
		stop("combinationsTest only works for binary data.")

	dd = rawData(obj)
	
	# select the variables we are interested in
	d = dd[,selected.vars]
	labels = dd$class
	
	if( length(selected.vars) == 1)
		comb = paste(d)
	else
		comb = apply(d, 1, paste, collapse="")
	
	# find out all combinations
	all.comb = (gtools::permutations(2, length(selected.vars), repeats.allowed=T)) - 1
	if( length(selected.vars) == 1){
		all.comb = paste(all.comb)
	} else {
		all.comb = apply(all.comb, 1, paste, collapse="")
	}

	counts = tapply(comb, labels, table)
	label.counts = table(labels)
	
	# do some fisher tests to see which differences are significant
	p.vals = rep(NA, length(all.comb))
	all.freq1 = rep(NA, length(all.comb))
	all.freq0 = rep(NA, length(all.comb)) 
	all.pct1 = rep(NA, length(all.comb))
	all.pct0 = rep(NA, length(all.comb))
	enrich.status = rep("", length(all.comb))
	enrich.fold = rep(0, length(all.comb))
	for(i in 1:length(all.comb)){
		cm = all.comb[i]
		if( !(cm %in% names(counts[["1"]]) ) )
			freq1 = c(0, label.counts[["1"]])
		else
			freq1 = c(counts[["1"]][[cm]], label.counts[["1"]] - counts[["1"]][[cm]])
			
		if( !(cm %in% names(counts[["0"]]) ) )
			freq0 = c(0, label.counts[["0"]])
		else
			freq0 = c(counts[["0"]][[cm]], label.counts[["0"]] - counts[["0"]][[cm]])
			
		p.vals[i] = fisher.test(rbind(freq1, freq0))$p.value
		all.freq1[i] = freq1[1]
		all.freq0[i] = freq0[1]
		all.pct1[i] = freq1[1] / (freq1[1]+freq1[2])
		all.pct0[i] = freq0[1] / (freq0[1]+freq0[2])
		if((freq1[1] / freq1[2]) > (freq0[1] / freq0[2])){
			enrich.status[i] = "enriched"
			enrich.fold[i] = round((freq1[1] / (freq1[1]+freq1[2])) / (freq0[1] / (freq0[1]+freq0[2]) ),2)
		} else if(freq1[1] == 0 & freq0[1] == 0){
			enrich.status[i] = NA
			enrich.fold[i] = NA
		} else {
			enrich.status[i] = "depleted"
			enrich.fold[i] = round((freq0[1] / (freq0[1]+freq0[2])) / (freq1[1] / (freq1[1]+freq1[2]) ),2)
		}
		if(verbose)
			cat(cm, p.vals[i], freq1[1], freq0[1], all.pct1[i], all.pct0[i], all.pct1[i]-all.pct0[i], enrich.status[i], enrich.fold[i], paste(which(comb==cm & labels==1), collapse=","), "\n")
	}
	
	p.vals = p.adjust(p.vals, method=p.adjust.method)
	
	# number of significant p vals
	sig = sum(p.vals < cutoff)

	# proportion of explained for 1 and 0
	exp1 = 0
	exp0 = 0	
	if(verbose)
		cat("\nMost important differences:\n")
	p.vals.order = order(p.vals)
	for(i in 1:length(p.vals.order)){
		j = p.vals.order[i]
		if(verbose){
			cat(all.comb[j], " ", p.vals[j], " ", all.freq1[j], all.freq0[j], all.pct1[j], all.pct0[j], all.pct1[j]-all.pct0[j], enrich.status[j], enrich.fold[j], "\n")
			if( i == sig )
				cat("----", cutoff, "significance cutoff\n")
		}
		if( i <= sig ){
			exp1 = exp1 + all.freq1[j]
			exp0 = exp0 + all.freq0[j]
		}
	}
	if(verbose){
		cat("\nExplained by significant differences:\n")
		cat("pos:", exp1, "out of", sum(labels), " ", exp1/label.counts[["1"]], "\n")
		cat("neg:", exp0, "out of", sum(!labels), " ", exp0/label.counts[["0"]], "\n")
		
		if(is.character(selected.vars))
			var.names.all = selected.vars
		else
			var.names.all = names(dd)[selected.vars]
		
		cat("\nJoint entropy of ", paste(var.names.all, collapse=", "), ": ", entropy(comb)*log2(exp(1)), " bits\n", sep="")
	}
	out = data.frame("combination"=all.comb, "p.value"=p.vals, "freq.pos"=all.freq1, "freq.neg"=all.freq0, 
		"type"=enrich.status, "fold.difference"=enrich.fold, stringsAsFactors=F)
	out = out[p.vals.order,]
	
	out
	
}


#' Calculate entropy from frequencies of observations for discrete data
#'
#' @param x the vector of frequencies, or a pdf of distribution
#' @return the entropy in bits
entropyFromFreq = function(x){
	x = x / sum(x) # make sure it is normalized to 1
	# calculate entropy in bits
	ent = sapply(x, function(y) ifelse(y==0, 0, -y*log2(y)))
	sum(ent)
}

#' Calculate the fold change when x is of size two (always show it >1)
#'
#' @param x input vector of size two
#' @return the proportion of x[1]/x[2] or x[2]/x[1] depending which is larger
foldChangeFromFreq = function(x){
	if( length(x) != 2 ){
		stop("foldChangeFromFreq requires a parameter of length two")
	}
	
	if( x[1] == 0 | x[2] == 0)
		return( NA )

	if( x[1] > x[2] )
		x[1] / x[2]
	else
		x[2] / x[1]
}


