\name{svmFeatureSelectionLOOCV}
\alias{svmFeatureSelectionLOOCV}
\title{Nested variable selection using LOOCV}
\usage{svmFeatureSelectionLOOCV(obj, selectionMode="direct", alpha=0.1, p.value.adjust.method="none",
    test.type="mc-x2", mc.replicates=5000, cost.range=logseq(0.01,
    1e+05, 8), gamma.range=logseq(1e-05, 100, 8), max.prop.SV=0.9,
    kernel="radial", skip.DDGraph=FALSE)
}
\description{Nested variable selection using LOOCV}
\details{A function to select variables in nested way using the following algorithm:
\enumerate{
\item repeat for each row in dataset:
\enumerate{
\item make new DDDataSet by removing one row and apply DDGraphs to select features
\item select best parameters using recalculateSVMparams (i.e. in an inner LOOCV loop)
\item make the classifier with best parameters and calculate output on the unseen row (removed in step 1)
}
\item return the collected predictions from step 1.3
}}
\value{the predictions for class labels from LOOCV}
\arguments{\item{obj}{the DDDataSet object}
\item{selectionMode}{which variables to take, possible values: 
"direct" (alias "p"), "direct and joint" (alias "ps"), "joint if no direct" (alias "snp")}
\item{alpha}{the alpha cutoff to use}
\item{p.value.adjust.method}{the p value adjustment for multiple testing to be applied}
\item{test.type}{the type of conditional independence test to be used}
\item{mc.replicates}{the number of Monte-Carlo replicates when determining p values}
\item{cost.range}{the range of cost parameter values to evaluate}
\item{gamma.range}{the range of gamma parameter values to evaluate}
\item{max.prop.SV}{the maximal proportion of support vectors to number of data points (rows in d)}
\item{kernel}{kernel type to use (takes valid package e1071 names like "radial")}
\item{skip.DDGraph}{if to skip DDGraph-based variable selection}
}

