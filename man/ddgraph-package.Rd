\name{ddgraph-package}
\alias{ddgraph-package}
\alias{ddgraph}
\docType{package}

\title{
  ddgraph package overview
}

\description{
    This package implements the Neighbourhood Consistent PC Algorithm (NCPC) for inferring the causal neighbourhood and Markov Blanket of a target variable, and a Direct Dependence Graphs (DDGraphs) for representing the conditional independence relationships. 

	The main goal of the NCPC algorithm is to infer direct from indirect dependencies of a set of variable to a target variable. The direct dependencies make up the causal 	neighbourhood of the target variable. This is achieved by performing conditional independence tests and therefore establishing statistical independence properties. 
	NCPC has been shown to have a larger recall rate in scenarios with highly correlated variables which are weakly associated to a sparse target variable. For more 
	details on the NCPC algorithm see (Stojnic et al, 2012). 
}

\details{
  \tabular{ll}{
    Package: \tab ddgraph\cr
    Type: \tab Package\cr
    License: \tab GPL-3\cr
    LazyLoad: \tab yes\cr
  }
  
  This package implements the NCPC/NCPC* algorithms, but also provides a unified front-end for inferring causal neighbourhood and Markov Blanket via Bayesian Network inference as provided by packages \code{bnlearn} and \code{pcalg}. 
  
  The package comes with two example datasets (Zizen et al 2009):
  \itemize{
  	\item \code{mesoBin} - binary dataset with 7 target variables - cis-regulatory module (CRM) classes. The variable correspond to transcription factor (TF) binding profiles over 1-5 time intervals. 
  	\item \code{mesoCont} - the original continuous version of the dataset. 
  }
  
  The main front-end function is \code{calcDependence()}. 
}

\author{

  Robert Stojnic\cr
  Cambridge Systems Biology Centre\cr
  University of Cambridge, UK

  Maintainer: Robert Stojnic \email{robert.stojnic@gmail.com}
}

\references{
	\itemize{
		\item R. Stojnic et al (2012): "A Graphical Modelling Approach to the Dissection of Highly Correlated Transcription Factor Binding Site Profiles", submitted.
		\item R. Zinzen et al (2009): "Combinatorial binding predicts spatio-temporal cis-regulatory activity"
Nature 462, no. 7269: 65-70.
	}
}

\keyword{ package }
