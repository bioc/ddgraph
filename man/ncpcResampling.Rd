\name{ncpcResampling}
\alias{ncpcResampling}
\title{NCPC Robustness from resampling}
\usage{ncpcResampling(obj, method="bootstrap", method.param, verbose=TRUE, ...)
}
\description{Estimate the NCPC robustness using either jackknife or bootstrap resampling.}
\details{Estimate the robustness of NCPC predictions (i.e. variable types: direct, joint, indirect, no dependence)
using resampling. Two type of resampling are available: bootstrap (where the whole dataset is resampled with 
replacement), and jackknifing (where 1 or more observation are removed at each resampling step). 

NCPC is run for the resampled datasets and statistics is produced about how many times is each variable
assigned one of the four types (direct, joint, indirect, no dependence). The final call for each variable is then
made according to the following algorithm (#direct is number of times variable is called direct):
\enumerate{
\item if #no dependence > #direct+joint+indirect => "no dependence"
\item else if #indirect > #direct+joint => "indirect"
\item else if #joint > #direct => "joint"
\item else "direct"
}}
\value{NCPCRobustness object with the raw results from resampling and summarized results}
\arguments{\item{obj}{the DDDataSet object}
\item{method}{the method to use to estimate how robust is the feature selection 
(valid values: "jackknife", or "bootstrap").}
\item{method.param}{the parameter to method, either number of data points to remove 
for "jackknife" (default: 1) or number of boostrap runs for "bootstrap" (default: 100).}
\item{verbose}{if to print out the progress}
\item{...}{other parameters to pass to ncpc()}
}
\examples{\dontrun{
# load the example data
data(mesoBin)

# run bootstrap resampling for NCPC with alpha=0.05
ncpcResampling(mesoBin$VM_SM, "bootstrap", 100, alpha=0.05)
# run bootstrap resampling for NCPC* with alpha=0.05
ncpcResampling(mesoBin$VM_SM, "bootstrap", 100, alpha=0.05, star=TRUE)

# run jackknifing for NCPC
ncpcResampling(mesoBin$VM_SM, "jackknife", 1, alpha=0.05)
}}

