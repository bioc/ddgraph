#
# Utility functions for generating random networks and difference testing scenarios
#


#' Estimate the in.degree distribution and conditional probability distribution from data
#'
#' The algorithm uses hill-climbing with BIC to construct the network and estimate the
#' parameters. Then, provided that for each in-degree there is at least two nodes, it
#' estimates the beta distribution parameters. 
#'
#' @title Estimate network distribution parameters
#' @param obj and object of class DDDataSet
#' @param use.class if to include the class variable into the estimate
#'
#' @return a list of two elements: in.degree.distr - distribution of in-degrees, and 
#'         beta.est - estimate beta distribution values
#' @export
#' @examples
#' data(mesoBin)
#' estimateNetworkDistribution(mesoBin$Meso)
estimateNetworkDistribution = function(obj, use.class=FALSE){

	tt = rawData(obj)
	if(!use.class)
		tt = tt[-which(names(tt)=="class")]
	
	orig.tt = tt
	
	if(dataType(obj) == "binary")
		tt = convertToFactor(tt)

	# generate the graph and the fitted graph with probabilities
	bn.graph = hc(tt)
	net = bn.fit(bn.graph, tt)

	# generate degree distribution
	all.degrees = rep(0, length(net))
	for(i in 1:length(net)){
		all.degrees[i] = length(dim(net[[i]]$prob))
	
	}

	# extract distributions of conditional probabilities for all degrees of freedom
	udegree = sort(unique(all.degrees))
	distr = list()

	for(i in 1:length(udegree)){
		inx = which(all.degrees == udegree[i])
		len = 2^udegree[i]
		m = matrix(0, ncol=len/2, nrow=length(inx))
		for(j in 1:length(inx)){
			m[j,] = as.vector(net[[ inx[j] ]]$prob)[seq(1,len,2)]
		}
	
		distr[[ udegree[i] ]] = m
	}

	# estimate parameters of the beta distribution
	beta.est = list()

	for(i in 1:length(distr)){
		if(i == 1){
			# use marginals for nodes when there is no parent
			val = colMeans(orig.tt)
			p = fitdistr(val, "beta", list(shape1=1, shape2=1))
			params = unlist(unlist(p)[c("estimate.shape1", "estimate.shape2")])
			params = data.frame("shape1" = params[1], "shape2" = params[2])
			beta.est[[i]] = params
		} else{
			# use the actual observed values, requires there is at least 2 values observed
			m = distr[[i]]
			params = matrix(0, ncol=2, nrow=ncol(m))
			for(j in 1:ncol(m)){
				val = m[,j]
				# values need to be in range (0,1) not exactly 1 or 0
				if(any(val == 1))
					val[val==1] = 0.99999
				if(any(val == 0))
					val[val==0] = 0.00001
			
				p = fitdistr(val, "beta", list(shape1=1, shape2=1))
				params[j,] = unlist(p)[c("estimate.shape1", "estimate.shape2")]
			}
			params = as.data.frame(params)
			names(params) = c("shape1", "shape2")
			beta.est[[i]] = params
		}
	}
	
	list(in.degree.distr = all.degrees-1, beta.est = beta.est)
}


#' Uniform distribution function for \code{random.bn.fit}
#'
#' Generate 2^n uniformly distributed numbers in range 0 to 1
#'
#' @title Uniform distribution for \code{random.bn.fit}
#' @param n number of variables
#' @return vector of 2^n uniform random numbers
#' @export
#' @examples
#' # return 8 random uniform numbers
#' prob.distr.unif(3)
prob.distr.unif = function(n){
	runif(2^n)
}

#' Generate 2^n numbers from distribution with most of the pdf mass in extreme probabilities (mirrored normal). 
#' We use standard deviation of 1/3 and modulo-1 of normal distribution. 
#'
#' @title Normal distribution function for \code{random.bn.fit}
#'
#' @param n number of variables
#' @param sd the standard deviation of distribution
#' @return vector of 2^n random numbers
#' @export
#' @examples
#' # return 8 random numbers since n=3
#' prob.distr.norm(3)
prob.distr.norm = function(n, sd=1/3){
	rnorm(2^n, sd=sd) %% 1
}

#' Convert graphNEL and friends representation to bn
#' 
#' @param graph graphNEL or graphAM object
graph.to.bn = function(graph){
	# start with an empty graph and add edges
	r = empty.graph(nodes(graph))
	
	p = as(graph, "graphAM")
	adj = p@adjMat
	rownames(adj) = colnames(adj)
	
	for(i in 1:nrow(adj)){
		for(j in 1:ncol(adj)){
			if(adj[i,j] == 1){
				r = set.arc(r, rownames(adj)[i], rownames(adj)[j])
			}
		}
	}
	
	r
	
}

#' Generate a random Bayesian network using package \code{bnlearn}. The nodes specify the partial ordering
#' of the graph, and the conditional probabilities are sampled from given distribution. The network is
#' generated to have on average given number of neighbours (i.e. both in-going and out-going edges)
#'
#' @title Generate a random \code{bn.fit} network
#'
#' @param nodes a vector of desired node names (basis for partial ordering)
#' @param num.neigh expected number of neighbours per node in the random graph
#' @param prob.distr the probability distribution function to use
#' @param bn.graph the \code{bn} object with an already laid out graph, if not supplied will be generated
#'
#' @return a list of two elements: \code{bn} - a \code{bn} object which contains the structure and 
#'         \code{bn.fit} - a \code{bn.fit} object with filled in conditional probabilities
#' @export
#' @examples
#' # a random network with 3 nodes "A", "B", "C" with average of 1 neighbour
#' random.bn.fit(c("A", "B", "C"), num.neigh=1)
random.bn.fit = function(nodes, num.neigh=2, prob.distr=prob.distr.norm, bn.graph=NULL){
	# generate a random graph structure using nodes for ordering
	if( is.null(bn.graph) ){
		#g = graph.to.bn(randomEGraph(nodes, num.neigh/(length(nodes)-1)))
		bn.graph = random.graph(nodes, prob=(num.neigh/(length(nodes)-1)))
	}

	# extract node information
	nodes = bn.graph$nodes	
	# start creating a bnfit object
	bnfit = lapply(nodes, function(x) list("parents"=x$parents, "children"=x$children))
	
	nodenames = names(nodes)
	for(i in 1:length(nodenames)){
		# set the name
		bnfit[[i]]$node = nodenames[i]
		
		n = length(bnfit[[i]]$parents)+1
		# generate some conditional probabilities
		prob = array(prob.distr(n), dim=rep(2,n))
		sel1 = c(1, rep("",n-1))
		sel2 = c(2, rep("",n-1))
		
		eval(parse(text=paste("prob[", paste(sel2, collapse=","), "] <- 1 - prob[", paste(sel1, collapse=","), "]", sep="")))
		
		# put in the names of dimensions
		dm = list()
		dm[[ nodenames[i] ]] = c("0", "1")
		if( n > 1 ){
			for(j in 1:(n-1)){
				dm[[ bnfit[[i]]$parents[j] ]] = c("0","1")
			}
		}
		
		dimnames(prob)  = dm
		class(prob) = "table"
		
		bnfit[[i]]$prob = prob
		
		class(bnfit[[i]]) = "bn.fit.dnode"
	}
	
	class(bnfit) = "bn.fit"
	
	list("bn"=bn.graph, "bn.fit"=bnfit)
		
}

#' Generate a random directed graph with the given node ordering and degree distribution
#'
#' @title Generate random network with degree distribution
#' @param nodes character vector of node names which species the node ordering
#' @param in.degree.distr the node in-degree distribution
#'
#' @return an object of class bn with the random graph
#' @export
#' @examples
#' # a random network of 5 nodes with choosen in-degree distribution
#' biased.graph(letters[1:5], c(0, 1, 1, 2, 2))
biased.graph = function(nodes, in.degree.distr){
	r = empty.graph(nodes)
	
	if( !any(in.degree.distr == 0) )
		stop("At least one node needs to have in-degree zero")
	
	stopifnot(length(nodes) == length(in.degree.distr))
	
	pool = in.degree.distr	
	for(i in 1:length(nodes)){
		num.parents = sample(pool[pool<i],1)
		# remove the picked value from the pool
		pool = pool[- (which(pool == num.parents)[1]) ] 
		if(num.parents == 0)
			next

		all.parents = sample(nodes[1:(i-1)], num.parents)		
		for(j in 1:num.parents){
			r = set.arc(r, all.parents[j], nodes[i])
		}
		
	}
	
	r
	
	
}

#' A version of random.bn.fit which generates a graph based on degree distribution and beta distribution for probabilities
#'
#' @title Random network with a biased degree distribution
#' @param nodes character vector of node names
#' @param beta.est the beta distribution parameters for different degrees of a node. Should be a list where [[2]] corresponds 
#'        to 2-dimenstional contingency table (i.e. one parent, one output). It contains a data.frame with columns shape1, shape2
#'        for the beta distribution, and rows are degrees of freedom (in this case 2, when P(Out=0|Parent=0) and P(Out=0|Parent=1))
#' @param in.degree.distr a vector with degree distribution for all the nodes in the network (names are ignored, and degree is 
#'        randomly sampled from this vector)
#' @param bn.graph if the graph structure is already available, then the graph structure in object of class "bn"  
#' @return a list of two elements: \code{bn} - a \code{bn} object which contains the structure and 
#'         \code{bn.fit} - a \code{bn.fit} object with filled in conditional probabilities
#' @export
#' @examples
#' # nodes, conditional probability distribution, an indegree distribution
#' nodes = letters[1:5]
#' beta.est = list(data.frame(shape1=2,shape2=3), data.frame(shape1=c(2,4), shape2=c(5,2)), data.frame(shape1=c(1,2,3,4), shape2=c(3,2,1,2)))
#' in.degree.distr = c(0, 1, 1, 2, 2)
#' # make a random graph using these parameters
#' biased.bn.fit(nodes, beta.est, in.degree.distr)
biased.bn.fit = function(nodes, beta.est, in.degree.distr, bn.graph=NULL){
	# generate a random graph structure using nodes for ordering
	if( is.null(bn.graph) ){
		bn.graph = biased.graph(nodes, in.degree.distr)
	}

	# extract node information
	nodes = bn.graph$nodes	
	# start creating a bnfit object
	bnfit = lapply(nodes, function(x) list("parents"=x$parents, "children"=x$children))
	
	nodenames = names(nodes)
	for(i in 1:length(nodenames)){
		# set the name
		bnfit[[i]]$node = nodenames[i]
		
		n = length(bnfit[[i]]$parents)+1
		# generate some conditional probabilities
		prob = rep(0, 2^(n-1))
		beta = beta.est[[n]]
		
		# there is 2^(n-1) degrees of freedom
		stopifnot(2^(n-1) == nrow(beta))		
		for(j in 1:nrow(beta)){
			prob[j] = rbeta(1, shape1=beta$shape1[j], shape2=beta$shape2[j])
		}
		# there is 2x more values that degrees of freedom
		prob = rep(prob, each=2)
		prob = array(prob, dim=rep(2,n))
		
		sel1 = c(1, rep("",n-1))
		sel2 = c(2, rep("",n-1))
		
		# make sure the probabilities are adding up to 1
		eval(parse(text=paste("prob[", paste(sel2, collapse=","), "] <- 1 - prob[", paste(sel1, collapse=","), "]", sep="")))
		
		# put in the names of dimensions
		dm = list()
		dm[[ nodenames[i] ]] = c("0", "1")
		if( n > 1 ){
			for(j in 1:(n-1)){
				dm[[ bnfit[[i]]$parents[j] ]] = c("0","1")
			}
		}
		
		dimnames(prob)  = dm
		class(prob) = "table"
		
		bnfit[[i]]$prob = prob
		
		class(bnfit[[i]]) = "bn.fit.dnode"
	}
	
	class(bnfit) = "bn.fit"
	
	list("bn"=bn.graph, "bn.fit"=bnfit)
		
}

#' Generate class labels by using the readout mechanism. Logical formula is applied to two variables
#' which are read out from the real data using the var1 and var2 probabilities. This only works
#' with binary variables. 
#'
#' @title Generate class labels by independent contributions of two variables
#'
#' @param data a matrix or data.frame containing binary observations (columns are variables)
#' @param var1 index or name of the first variable
#' @param var2 index or name of the second variable
#' @param var1.prob1 the conditional probability P(class labels = 1|var1=1)
#' @param var1.prob0 the conditional probability P(class labels = 1|var1=0)
#' @param var2.prob1 the conditional probability P(class labels = 1|var2=1)
#' @param var2.prob0 the conditional probability P(class labels = 1|var2=0)
#' @param logical.formula logical formula to apply
#' @param false.neg a false negative probability
#' @param false.pos a false positive probability
#' @return a binary vector containing the class labels
#' @export
#' @examples
#' # noisy OR function with 0.1 probability of error for reading "a" and "b" (error in both 1 and 0)
#' data <- cbind("a"=c(0,0,1,1), "b"=c(0,1,0,1))
#' independent.contributions.formula(data, "a", "b", 0.9, 0.1, 0.9, 0.1, "a | b")
independent.contributions.formula = function(data, var1, var2, var1.prob1, var1.prob0, var2.prob1, var2.prob0, logical.formula,
	false.neg=0, false.pos=0){

	apply.formula = function(x){
		expr = paste("with(as.list(x),", logical.formula, ")", sep="")
		out = eval(parse(text=expr))
		as.numeric(out)
	}
	
	# extract the two variables
	a = data[,var1]
	b = data[,var2]
	
	labels = rep(0, nrow(data))
	
	for(i in 1:nrow(data)){
		out = c(0,0)
		names(out) = colnames(data[,c(var1,var2)])
		
		if( a[i] == 1){
			if(runif(1) <= var1.prob1)
				out[1] = 1
		}
		if( a[i] == 0){
			if(runif(1) <= var1.prob0)
				out[1] = 1
		}

		if( b[i] == 1){
			if(runif(1) <= var2.prob1)
				out[2] = 1
		}
		if( b[i] == 0){
			if(runif(1) <= var2.prob0)
				out[2] = 1
		}

		labels[i] = apply.formula(out)
		
		# apply false negative rate
		if(labels[i] == 1){
			if(runif(1) < false.neg)
				labels[i] = 0
		}
		# false positive rate
		if(labels[i] == 0){
			if(runif(1) < false.pos)
				labels[i] = 1
		}
	
	}

	labels	
}

#' Generate class labels by using the readout mechanism. Logical formula is applied to two variables
#' which are read out from the real data using the var1 and var2 probabilities. This only works
#' with binary variables. 
#'
#' @title Generate class labels by a noisy formula with high false negative rate
#'
#' @param data a matrix or data.frame containing binary observations (columns are variables)
#' @param var1 index or name of the first variable
#' @param var2 index or name of the second variable
#' @param false.neg a false negative probability
#' @param logical.formula logical formula to apply
#' @return a binary vector containing the class labels
#' @export
#' @examples
#' # noisy OR function with 0.1 probability of error for reading "a" and "b" (error in both 1 and 0)
#' data <- cbind("a"=c(0,0,1,1), "b"=c(0,1,0,1))
#' formulaFalseNeg(data, "a", "b", 0.8, "a | b")
formulaFalseNeg = function(data, var1, var2, false.neg, logical.formula){

	apply.formula = function(x){
		expr = paste("with(as.list(x),", logical.formula, ")", sep="")
		out = eval(parse(text=expr))
		as.numeric(out)
	}
	
	# extract the two variables
	a = data[,var1]
	b = data[,var2]
	
	labels = rep(0, nrow(data))
	
	for(i in 1:nrow(data)){
		out = c(a[i], b[i])
		names(out) = colnames(data[,c(var1,var2)])
		
		labels[i] = apply.formula(out)
		
		# apply false negative rate
		if(labels[i] == 1){
			if(runif(1) < false.neg)
				labels[i] = 0
		}
	
	}

	labels	
}


#' Version of \code{independent.contributions.formula} that works with any number of variables. See the help page
#' for \code{independent.contributions.formula} for description of functionality. 
#'
#' @title Generate class labels by independent contributions of two variables
#'
#' @param data a matrix or data.frame containing binary observations (columns are variables)
#' @param target.vars indexes of target variables
#' @param prob1 vector of P(class labels = 1|varX=1) for different X
#' @param prob0 vector of P(class labels = 1|varX=0) for different X
#' @param logical.formula a character string for the formula
#'
#' @return a vector of binary class labels
#' @export
#' @examples
#' # noisy OR function with three variables and with noise level of 0.1 for a, b, and 0.2 for c
#' data <- cbind("a"=c(0,0,0,0,1,1,1,1), "b"=c(0,0,1,1,0,0,1,1), "c"=c(0,1,0,1,0,1,0,1))
#' independent.contributions.formula.mul(data, c("a", "b", "c"), c(0.9, 0.9, 0.8), c(0.1, 0.1, 0.2), "a | b | c")
independent.contributions.formula.mul = function(data, target.vars, prob1, prob0, logical.formula){
	apply.formula = function(x){
		expr = paste("with(as.list(x),", logical.formula, ")", sep="")
		out = eval(parse(text=expr))
		as.numeric(out)
	}
	
	# extract the two variables
	vars = data[,target.vars]
	
	labels = rep(0, nrow(data))
	
	for(i in 1:nrow(data)){
		out = rep(0, ncol(vars))
		names(out) = colnames(data[,target.vars])
		
		for(j in 1:ncol(vars)){
			# probability of outputing 0 if input is one
			if(vars[i,j] == 1){
				if(runif(1) <= prob1[j])
					out[j] = 1
			}
			# probability of outputing 1 if input is zero
			if(vars[i,j] == 0){
				if(runif(1) <= prob0[j])
					out[j] = 1
			}
		}

		labels[i] = apply.formula(out)
	
	}

	labels	
}


