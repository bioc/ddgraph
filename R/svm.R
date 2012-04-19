#
# Support Vector Machine (SVM) functions for asssessing the performance (and selecting variables) using SVMs 
# 


#' Generate sequence but in log scale. This function takes takes the length of log-sequence and the minimal and
#' maximal point. It returns the interval between \code{a} and \code{b} divided in log scale. 
#'
#' @title Generate sequence in log scale
#' 
#' @param a the smaller value in the interval
#' @param b the bigger value in the interval
#' @param n the number of intervals to divide a,b into 
#'
#' @return a vector of numbers
#' @export
#' @examples
#' # produces vector c(0.01, 0.1, 1)
#' logseq(0.01, 1, 3)
logseq = function(a,b,n=8) exp(seq(log(a), log(b), length.out=n))


#' Nested variable selection using LOOCV
#'
#' A function to select variables in nested way using the following algorithm:
#' \enumerate{
#'   \item repeat for each row in dataset:
#'   \enumerate{
#'      \item make new DDDataSet by removing one row and apply DDGraphs to select features
#'      \item select best parameters using recalculateSVMparams (i.e. in an inner LOOCV loop)
#'      \item make the classifier with best parameters and calculate output on the unseen row (removed in step 1)
#'   }
#'   \item return the collected predictions from step 1.3
#' }
#'
#' @title Nested variable selection using LOOCV
#'
#' @param obj the DDDataSet object
#' @param selectionMode which variables to take, possible values: 
#'        "direct" (alias "p"), "direct and joint" (alias "ps"), "joint if no direct" (alias "snp")
#' @param alpha the alpha cutoff to use
#' @param p.value.adjust.method the p value adjustment for multiple testing to be applied
#' @param test.type the type of conditional independence test to be used
#' @param mc.replicates the number of Monte-Carlo replicates when determining p values
#' @param cost.range the range of cost parameter values to evaluate
#' @param gamma.range the range of gamma parameter values to evaluate
#' @param max.prop.SV the maximal proportion of support vectors to number of data points (rows in d)
#' @param kernel kernel type to use (takes valid package e1071 names like "radial")
#' @param skip.DDGraph if to skip DDGraph-based variable selection
#'
#' @return the predictions for class labels from LOOCV
svmFeatureSelectionLOOCV = function(obj, selectionMode="direct", alpha=0.1, p.value.adjust.method="none", 
    test.type="mc-x2", mc.replicates=5000, cost.range = logseq(1e-2, 1e+5, 8), gamma.range = logseq(1e-5, 1e2, 8),
    max.prop.SV = 0.9, kernel = "radial", skip.DDGraph=FALSE){
	
	# algorithm outline:
	# 1) repeat for each row in dataset:
	#   1.1) make new DDDataSet by removing one row and apply DDGraphs to select features
	#   1.2) select best parameters using recalculateSVMparams (i.e. in an inner LOOCV loop)
	#   1.3) make the classifier with best parameters and calculate output on the unseen row (removed in step 1)
	# 2) return the collected predictions from step 1.3
	
	d = rawData(obj)
	d.fac = convertToFactor(d)
	
	# the predicted class values
	predicted.value = rep(NA, nrow(d))
	selected.vars = list()
	selected.params = list()
	
	# cache values as many data frames in LOOCV can be the same
	cache = list()
	
	for(i in 1:nrow(d)){
		message("\nRemoving row ", i, " / ", nrow(d))
		# make new dataset by deleting a row
		nd = makeDDDataSet(d[-i,], name=paste(datasetName(obj), "without row", i))
		
		if(!skip.DDGraph){		
			# make cdgraphs
			eg = ncpc(nd, alpha=alpha, p.value.adjust.method=p.value.adjust.method, test.type=test.type, mc.replicates=mc.replicates)
		
			# select variables
			if(selectionMode %in% c("direct", "p")){
				vars = eg$direct
			} else if(selectionMode %in% c("direct and joint", "ps")){
				vars = eg$directAndJoint
			} else if(selectionMode %in% c("joint if no direct", "snp")){
				vars = eg$direct
				if(length(vars)==0)
					vars = eg$joint
			}
		} else{
			# take all variables
			vars = 1:(ncol(d)-1)
		}
		
		class.weight = 1/table(convertToFactor(nd$class))
		
		# select the variables
		seldf = data.frame(rawData(nd)[,vars], nd$class)
		names(seldf) = c(names(nd)[vars], "class")
		
		message("\nSelecting cost and gamma parameters based on inner LOOCV")
		message("Using ", length(vars) ," variable(s): ", paste(names(nd)[vars], sep=", "), "\n")
		# first try to get the parameters from cache
		from.cache = FALSE
		if( length(cache) > 0){
			for(j in 1:length(cache)){
				# match the matrices, if they are the same, the results are going to be the same as well
				if( all(dim(cache[[j]]$seldf) == dim(seldf)) ){
					if( all(cache[[j]]$seldf == seldf) ){
						params = cache[[j]]$params
						from.cache = TRUE
						message("Getting parameters from cache...")
						break
					}
				}
			}
		}
		if(!from.cache){
			# select parameters based on these variables
			params = recalculateSVMparams(cost.range, gamma.range, convertToFactor(seldf), class.weight, kernel, max.prop.SV)
			cache[[length(cache)+1]] = list(seldf = seldf, params = params)
		}
		
		message("Using parameters cost = ", params$cost, ", gamma = ", params$gamma, " with AUC = ", params$auc, "\n")
		
		# construct the model
		model = svm(class~., data=convertToFactor(seldf), cost=params$cost, gamma=params$gamma, class.weight=class.weight, kernel=kernel, type="C-classification")
				
		# write out the results
		predicted.value[i] = predSVM(model, d.fac[i,c(vars, which(names(d.fac)=="class"))])
		selected.vars[[i]] = vars
		selected.params[[i]] = params
	}
	
	list("predicted.value"=predicted.value, "vars"=selected.vars, "params"=selected.params)
}

#' A companion function for \code{svmFeatureSelectionLOOCV()} to plot the results. 
#'
#' @title Plot SVM performance into a pdf file
#'
#' @param obj the DDDataSet object for which the SVM performance is measured
#' @param results the results from svmFeatureSelectionLOOCV
#' @param plot.file the name of the output pdf file
plotSVMPerformance = function(obj, results, plot.file=NULL){
	# extract results
	predicted.class = results$predicted.value
	vars = results$vars
	target.class = obj$class
	
	# use ROCR to get the metrics
	pred = prediction(round(predicted.class,8), ifelse(target.class==0,-1,1))
	# print out different metrics and save to plots
	auc = performance(pred, "auc")@y.values[[1]]
	message("AUC: ", auc)
	
	# get the variable statistics
	vars.stat = sort(table(unlist(vars)))
	names(vars.stat) = names(obj)[as.numeric(names(vars.stat))]
	vars.stat = vars.stat / nrow(rawData(obj))	
	
	# get stats on the parameter usage
	gamma.stat = sort(table(sapply(results$params, function(x) x$gamma))) / nrow(rawData(obj))
	cost.stat = sort(table(sapply(results$params, function(x) x$cost))) / nrow(rawData(obj))
	auc.stat = sapply(results$params, function(x) x$auc)
	
	if( !is.null(plot.file) )
		pdf(plot.file, height=7*3, width=7*2, pointsize=22)
	par(mfrow=c(3,2))
	# PRC plot
	plot(performance(pred,"prec","rec"), col="blue", xlim=c(0,1), ylim=c(0,1), main="PRC plot from LOOCV")
	# ROC
	plot(performance(pred,"tpr","fpr"), col="blue", main=paste("ROC with AUC =", round(auc,3)))
	abline(0,1)
	
	barplot(vars.stat, horiz=TRUE, main="Proportion of each variable being choosen", xlab="proportion of total LOOCV runs")
	
	barplot(gamma.stat, horiz=TRUE, main="Gamma parameter", xlab="proportion of total LOOCV runs")
	barplot(cost.stat, horiz=TRUE, main="Cost parameter", xlab="proportion of total LOOCV runs")	
	hist(auc.stat, breaks=20, main="AUC distribution", freq=FALSE, xlab="AUC from LOOCV")
	
	if( !is.null(plot.file) )
		dev.off()
}

#' Leave-one-out cross validation systematically leaves out one row from the data, retrains the
#' classifier and then uses the retrained classifier to make a prediction for the left-out row.
#'
#' @title Leave-one-out cross validation
#'
#' @param data The \code{data.frame} with data. Columns are variables, rows are observations.
#' @param train.fun The training function that takes the data without one of the rows left out.
#' @param eval.fun The prediction function that takes the trained model and the left out data point.
#' @param verbose If to print progress indication
#'
#' @return A vector of length \code{nrow(data)} containing predictions from \code{eval.fun} when each 
#'        row is left out once
loocv = function(data, train.fun, eval.fun, verbose=FALSE){	
	score = rep(0, nrow(data))
	
	for(i in 1:nrow(data)){
		if(verbose)
			message(".", appendLF=FALSE)
		subset = setdiff( 1:nrow(data), i )
		stopifnot( is.factor(data$class) )
		# apply the fitting function on subset of data
		model = train.fun(data[subset,])
		score[i] = eval.fun(model, data[i,])
	}
	if(verbose)
		message("")
	
	score
	
}

#' Calculate the decision value of an SVM model. Note this is different from the actual
#' prediction which is either 0 or 1, while decision values go from -1 to 1.
#' (taken from [Zizen 2009] supplementary code)
#'
#' @title Calculate the decision value of an SVM model
#'
#' @param f The trained SVM model object.
#' @param feature The input value to which output is needed.
#'
#' @return Decision value in the range -1 to 1.
predSVM = function(f, feature){
    as.vector( - diff(f$labels)*attr( predict(f, feature, decision.values=TRUE ), 
	"decision.values") )
}


#' Find the cost/gamma parameters based on a grid search by best AUC and by 
#' limiting the number of support vectors. Currently only supports discreet binary data. 
#'
#' @title Calculate SVM hyperparameters based on grid search
#'
#' @param cost.range the range of cost parameter values to evaluate
#' @param gamma.range the range of gamma parameter values to evaluate
#' @param d the data.frame with variables as columns, the class labels must be labelled with "class"
#' @param class.weight the class weights to use (if there is an large bias for positive/negative class)
#' @param kernel kernel type to use (takes valid package e1071 names like "radial")
#' @param max.prop.SV the maximal proportion of support vectors to number of data points (rows in d)
#' 
#' @return a list with the two parameters that give best AUC in LOOCV
#' @export
#' @examples
#' \dontrun{
#' data(mesoBin)
#' # get SVM AUC etc over cost rage of 1, 100, and gamma range of 0.1, 1
#' recalculateSVMparams(c(1, 100), c(0.1, 1), convertToFactor(rawData(mesoBin$Meso)))
#' }
recalculateSVMparams = function(cost.range, gamma.range, d, 
                                class.weight=1/table(convertToFactor(d$class)), 
                                kernel="radial", max.prop.SV=0.9){
	svm.params = expand.grid(cost=cost.range, gamma=gamma.range)
	svm.params.auc = rep(0, nrow(svm.params))
	svm.params.normal.nSV = rep(T, nrow(svm.params))
	# record AUC for various parameter combinations
	for(iteration in 1:2){
		if(iteration == 2)
			message("Entering second iteration, relaxing maximal SV number constrain")
		for(i in 1:nrow(svm.params)){
			message("Trying cost = ", svm.params$cost[i], ", gamma = ", svm.params$gamma[i])
		
			stopifnot( is.factor(d$class) )
		
			# check the number of support vectors
			svm.fit = svm(class~., data=d, cost=svm.params$cost[i], gamma=svm.params$gamma[i], 
                          class.weight=class.weight, kernel=kernel, type="C-classification")			
			nSV = sum(svm.fit$nSV)
			# in first iteraction check if SV number is not too large
			if( (nSV / nrow(d)) > max.prop.SV & iteration == 1){
				message("Too many support vectors, skipping...")
				svm.params.normal.nSV[i] = F
				next
			}
		
			ml = loocv(d, function(dd){ svm(class~., data=dd, cost=svm.params$cost[i], type="C-classification", 
				gamma=svm.params$gamma[i], class.weight=class.weight, kernel=kernel) }, function(m,d){ predSVM(m, d) })
			
			# take only the first 4 significant decimals, as taking all of them might introduce numerical artefacts!
			ml = round(ml,4)
				
			pred.ml = prediction(ml, ifelse(d$class==0, -1, 1))	
			svm.params.auc[i] = performance(pred.ml, "auc")@y.values[[1]]
			message("AUC = ", svm.params.auc[i], "\n")
		}
	
		# if we have one one which is not over-fitting then we are finished
		if( iteration == 1 & any(svm.params.normal.nSV) )
			break
	}
	
	# select the param with max AUC
	best.inx = which.max(svm.params.auc)
	
	list( "cost" = svm.params$cost[best.inx], "gamma" = svm.params$gamma[best.inx], 
    "auc"=svm.params.auc[best.inx], "svm.params"=cbind(svm.params,
                                                       "auc"=svm.params.auc,"nSV"=svm.params.normal.nSV) )
}


