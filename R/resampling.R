#
# Resampling using jackknife and bootstrap for NCPC
#


#' Estimate the NCPC robustness using either jackknife or bootstrap resampling. 
#'
#' Estimate the robustness of NCPC predictions (i.e. variable types: direct, joint, indirect, no dependence)
#' using resampling. Two type of resampling are available: bootstrap (where the whole dataset is resampled with 
#' replacement), and jackknifing (where 1 or more observation are removed at each resampling step). 
#'
#' NCPC is run for the resampled datasets and statistics is produced about how many times is each variable
#' assigned one of the four types (direct, joint, indirect, no dependence). The final call for each variable is then
#' made according to the following algorithm (#direct is number of times variable is called direct):
#' \enumerate{
#'   \item if #no dependence > #direct+joint+indirect => "no dependence"
#'   \item else if #indirect > #direct+joint => "indirect"
#'   \item else if #joint > #direct => "joint"
#'   \item else "direct"
#' }
#'
#' @title NCPC Robustness from resampling
#' @param obj the DDDataSet object
#' @param method the method to use to estimate how robust is the feature selection 
#'        (valid values: "jackknife", or "bootstrap").
#' @param method.param the parameter to method, either number of data points to remove 
#'        for "jackknife" (default: 1) or number of boostrap runs for "bootstrap" (default: 100).
#' @param verbose if to print out the progress
#' @param ... other parameters to pass to ncpc()
#'
#' @return NCPCRobustness object with the raw results from resampling and summarized results
#' @export
#' @examples
#' \dontrun{
#' # load the example data
#' data(mesoBin)
#' 
#' # run bootstrap resampling for NCPC with alpha=0.05
#' ncpcResampling(mesoBin$VM_SM, "bootstrap", 100, alpha=0.05)
#' # run bootstrap resampling for NCPC* with alpha=0.05
#' ncpcResampling(mesoBin$VM_SM, "bootstrap", 100, alpha=0.05, star=TRUE)
#' 
#' # run jackknifing for NCPC
#' ncpcResampling(mesoBin$VM_SM, "jackknife", 1, alpha=0.05)
#' }
ncpcResampling = function(obj, method="bootstrap", method.param=NULL, verbose=TRUE, ...){
	# consistency checking
	if( !(method %in% c("bootstrap", "jackknife")) ){
		stop(paste("Unknown method \"", method, "\" selected, valid values are: jackknife, bootstrap", sep=""))
	}
	
	if(!is.null(method.param) && !is.numeric(method.param))
		stop("Error: parameter \"method.param\" needs to be either NULL or a numeric value")
	
	# record the input parameters for this function
	params = list()
	for(p in names(formals())){
		if( p == "" | p == "..." | p == "obj")
			next
		
		params[[p]] = eval(parse(text=p))
	}
	
	# raw dataframe with the variables
	d = rawData(obj)
	nm = names(d)
	
	# find out the method
	if(method == "bootstrap"){
		if( is.null(method.param) )
			method.param = 100
		runs = method.param
		# make the selection matrix for bootstrap
		sel.matrix = matrix(sample(nrow(d), nrow(d)*runs, replace=T), nrow=runs)
		
	} else if(method == "jackknife"){
		if( is.null(method.param) )
			method.param = 1
		if( nrow(d) %% method.param != 0)
			stop("The number of jackknife data points to remove need to be an exact multiple to number of data points")
		runs = nrow(d) / method.param
		# the selection matrix for jackknife
		sel.matrix = - matrix(sample(nrow(d), nrow(d), replace=F), ncol=method.param)
		
	}	
	
	## the results list
	res = list()
	
	# do the jackknife/boostrap runs
	for(i in 1:runs){
		if(verbose)
			message("\n\nRun ", i, " / ", runs)
		
		sel = sel.matrix[i,]			
		nd = makeDDDataSet(d[sel,], name=paste(datasetName(obj), method, "run", i))
		
		# consistency checking, expect certain number of rows of data..
		if(method == "bootstrap")
			stopifnot(nrow(rawData(nd)) == nrow(d))
		if(method == "jackknife")
			stopifnot(nrow(rawData(nd)) == (nrow(d)-method.param))
			
		# make DDGraphs!
		eg = ncpc(nd, verbose=verbose, ...)
		
		if(i == 1){
			# append the parameter sets
			params = c(params, eg$params)
		}
		
		# extract the direct/joint/indirect etc		
		if( is.null(eg) )
			res[[i]] = list("direct"=NULL, "joint"=NULL, "indirect"=NULL, "conditional"=NULL,
				"ps"=NULL, "snp"=NULL)
		else{
			if( length(eg$direct) == 0 )
				snp = eg$joint
			else
				snp = eg$direct
				
			res[[i]] = list("direct"=nm[eg$direct], "joint"=nm[eg$joint], "indirect"=nm[eg$indirect],
				"conditional"=nm[eg$conditional], "directAndJoint"=nm[eg$directAndJoint], "jointIfNotDirect"=nm[snp])
				
			if(dataType(nd) == "binary")
				res[[i]][["positiveClassSize"]] = sum(nd$class)
		}
				
	}
	
	# make the resulting object
	calculateNCPCRobustnessStats(  makeNCPCRobustness(obj, res, params) )
	
}

#' Make a new NCPCRobustness object 
#'
#' Make a new NCPCRobustness object just with the raw resampling data and parameters used to generate them.
#' Should never directly use this function, but only via \code{DDDataSet::NCPCRobustness()}.
#'
#' @param dataset the DDDataSet object
#' @param raw the list of raw resampling classification of variables (direct, joint, etc..)
#' @param params the parameters used to generate the data (only the non-default one are listed)
#'
#' @return a new NCPCRobustness object
makeNCPCRobustness = function(dataset, raw, params){
	new("NCPCRobustness", dataset=dataset, raw=raw, params=params, tables=list(), runs=length(raw),
		enriched.pss = data.frame(), enriched.ps=data.frame(), not.enriched=data.frame(), final.calls=data.frame())	
}

#' Calculate NCPCRobustness statistics 
#'
#' Calculate the statistics for the NCPCRobustness object - this is separate from object 
#' construction for convenience of testing, should always be called after object creation. 
#' Never use directly (except for testing), use instead via \code{DDDataSet::NCPCRobustness()}.
#'
#' @param obj NCPCRobustness object
#' @return the modified NCPCRobustness object with the statistics calculated
calculateNCPCRobustnessStats = function(obj){	
	# get the raw results
	res = obj@raw
	
	# make the sorted tables
	for(name in c("direct", "joint", "indirect", "directAndJoint", "jointIfNotDirect", "conditional")){
		obj@tables[[name]] = sort(table(unlist(lapply(res, function(x) x[[name]]))), decreasing=TRUE)
	}
	obj@tables[["positiveClassSize"]] = table(unlist(lapply(res, function(x) x[["positiveClassSize"]])))
	
	var.names = variableNames(obj@dataset)
	
	# extract the counts in consistent ordering
	direct = rep(0, length(var.names))
	indirect = rep(0, length(var.names))
	joint = rep(0, length(var.names))
	conditional = rep(0, length(var.names))
	
	direct[match( names(obj@tables$direct), var.names)] = obj@tables$direct / obj@runs
	joint[match( names(obj@tables$joint), var.names)] = obj@tables$joint / obj@runs
	indirect[match( names(obj@tables$indirect), var.names)] = obj@tables$indirect / obj@runs
	conditional[match( names(obj@tables$conditional), var.names)] = obj@tables$conditional / obj@runs		
	not.enriched = round(1 - direct - indirect - joint - conditional, 10) # weed out any numerical error with round()

	# split variables into enriched or not	
	freq.enriched = cbind(direct+joint+indirect+conditional, not.enriched)
	
	enrich.calls = apply(freq.enriched, 1, which.max)
	enrich.entropy = apply(freq.enriched, 1, entropyFromFreq)
	enrich.fc = apply(freq.enriched, 1, foldChangeFromFreq)	
	# selectors for enriched and non-enriched variables
	sel.enrich = (enrich.calls == 1)
	sel.not.enrich = (enrich.calls == 2)
	
	stopifnot( all(sel.enrich | sel.not.enrich) ) # consistency check
	
	# put together all the data and select the no dependence ones
	enrich.all = data.frame(freq.enriched, enrich.fc, enrich.entropy)
	names(enrich.all) = c("informative", "no dependence", "fold difference", "uncertainty")
	rownames(enrich.all) = var.names
	obj@not.enriched = enrich.all[sel.not.enrich,]
	
	# split the enriched ones by direct/indirect/joint
	freq.pss = cbind(direct, joint, indirect)
	pss.entropy = apply(freq.pss, 1, entropyFromFreq)
	pss.calls = c("direct", "joint", "indirect")[apply(freq.pss, 1, which.max)]
	
	freq.ps = cbind(direct+joint, indirect)
	ps.entropy = apply(freq.ps, 1, entropyFromFreq) 
	ps.fc = apply(freq.ps, 1, foldChangeFromFreq)
	ps.calls = c("direct and joint", "indirect")[apply(freq.ps, 1, which.max)]
	
	enriched.pss = data.frame(freq.pss, pss.calls, pss.entropy, enrich.fc, enrich.entropy)
	enriched.ps = data.frame(freq.ps, ps.calls, ps.fc, ps.entropy, enrich.fc, enrich.entropy)
	
	names(enriched.pss) = c("direct", "joint", "indirect", "final call", "p/s/s uncertainty", 
		"enrichment fold difference", "enrichment uncertainty")
	names(enriched.ps) = c("direct and joint", "indirect", "final call", "ps/s fold difference", 
		"ps/s uncertainty", "enrichment fold difference", "enrichment uncertainty")
	rownames(enriched.pss) = var.names
	rownames(enriched.ps) = var.names
	
	obj@enriched.pss = enriched.pss[sel.enrich,]
	obj@enriched.ps = enriched.ps[sel.enrich,]
	
	# order by entropy
	obj@not.enriched = obj@not.enriched[ order(obj@not.enriched[,4]), ]
	obj@enriched.pss = obj@enriched.pss[ order(obj@enriched.pss[,5]), ]
	obj@enriched.ps = obj@enriched.ps[ order(obj@enriched.ps[,5]), ]	
	
	# make the final calls according to following algorithm:
	# 1) if not.enriched > enriched => not.enriched
	# 2) else if indirect > direct+joint => indirect
	# 3) else if joint > direct => joint
	# 4) else direct
	
	final.calls = data.frame("name"=var.names, "type"="", "probability"=0, stringsAsFactors=F)
	
	for(i in 1:nrow(final.calls)){
		if( not.enriched[i] > (1-not.enriched[i]) ){
			final.calls$type[i] = "no dependence"
			final.calls$probability[i] = not.enriched[i]
		} else if( conditional[i] > direct[i] + indirect[i] + joint[i] ){
			final.calls$type[i] = "conditional"
			final.calls$probability[i] = conditional[i]
		} else if( indirect[i] > direct[i] + joint[i] ){
			final.calls$type[i] = "indirect"
			final.calls$probability[i] = indirect[i]
		} else if( joint[i] > direct[i] ){
			final.calls$type[i] = "joint"
			final.calls$probability[i] = joint[i]
		} else {
			final.calls$type[i] = "direct"
			final.calls$probability[i] = direct[i]
		}
	}
	
	obj@final.calls = final.calls
		
	obj
}

