\name{NCPCRobustness-class}
\alias{NCPCRobustness-class}
\alias{NCPCRobustness}
\title{NCPC resampling robustness...}
\description{NCPC resampling robustness}
\details{Data class that stores the robustness information associated with an NCPC result from
resampling runs (bootstrap of jackknifing). It contains the results from the resampling
runs as well as summary statistics. The \code{final.calls} slot contains the final
assigned types based on resampling.}
\section{Slots}{\describe{\item{\code{dataset}:}{(\code{\link{DDDataSet}}) the associated DDDataSet object}\item{\code{raw}:}{(\code{\link{list}}) the raw data from the robustness analysis}\item{\code{params}:}{(\code{\link{list}}) the parameters used to generate the data (including the resampling method)}\item{\code{tables}:}{(\code{\link{list}}) the frequencies of assigning each variable to a class}\item{\code{runs}:}{(\code{\link{numeric}}) the number of resampling runs}\item{\code{enriched.pss}:}{(\code{\link{data.frame}}) the table with reports for consistently enriched variables split}\item{\code{enriched.ps}:}{(\code{\link{data.frame}}) the table with reports for consistently enriched variable split
into two classes: directAndJoint, indirect}\item{\code{not.enriched}:}{(\code{\link{data.frame}}) the table with reports for the consistently not enriched variables}\item{\code{final.calls}:}{(\code{\link{data.frame}}) the table with finals calls for types of variables}}}

