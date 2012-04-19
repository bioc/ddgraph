\name{recalculateSVMparams}
\alias{recalculateSVMparams}
\title{Calculate SVM hyperparameters based on grid search}
\usage{recalculateSVMparams(cost.range, gamma.range, d,
    class.weight=1/table(convertToFactor(d$class)), kernel="radial",
    max.prop.SV=0.9)
}
\description{Find the cost/gamma parameters based on a grid search by best AUC and by 
limiting the number of support vectors. Currently only supports discreet binary data.}
\value{a list with the two parameters that give best AUC in LOOCV}
\arguments{\item{cost.range}{the range of cost parameter values to evaluate}
\item{gamma.range}{the range of gamma parameter values to evaluate}
\item{d}{the data.frame with variables as columns, the class labels must be labelled with "class"}
\item{class.weight}{the class weights to use (if there is an large bias for positive/negative class)}
\item{kernel}{kernel type to use (takes valid package e1071 names like "radial")}
\item{max.prop.SV}{the maximal proportion of support vectors to number of data points (rows in d)}
}
\examples{\dontrun{
data(mesoBin)
# get SVM AUC etc over cost rage of 1, 100, and gamma range of 0.1, 1
recalculateSVMparams(c(1, 100), c(0.1, 1), convertToFactor(rawData(mesoBin$Meso)))
}}

