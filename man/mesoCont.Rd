\name{mesoCont}
\alias{mesoCont}
\docType{data}

\title{
  A list of continuous \code{DDDataSet} objects.
}
\description{
  \code{mesoCont} is a list of objects of class \code{\linkS4class{DDDataSet}}. 
  It has been generated with the following code:
  
  \code{mesoCont <- toDDDataSet(readFurlongData(), prettyNames=TRUE, convertToBinary=FALSE)}
}

\usage{data(mesoCont)}

\details{
The dataset represents original continuous binding signal for 5 transcription factors (TFs) at 1-5 time points
during embryonic mesoderm development in Drosophila Melanogaster (Zinzen et al, 2009). The
original data is retained (from Supplementary Table 8 of the paper).  

The list contains 7 objects of type \code{\linkS4class{DDDataSet}} for 7 cis-regulatory
module (CRM) classes. These classes are: neg (negative class of CRMs), Meso (CRMs active in
early mesoderm), Meso_SM (CRMs active in early mesoderm and somatic muscle), VM (visceral 
muscle), SM (somatic muscle), VM_SM (active in both somatic and visceral muscle) and CM
(active in cardiac muscle). 
}

\examples{
data(mesoCont)
names(mesoCont)
class(mesoCont$VM)
}
\references{Robert P. Zinzen et al., "Combinatorial binding predicts spatio-temporal cis-regulatory activity,"
Nature 462, no. 7269 (November 5, 2009): 65-70.}
\keyword{datasets}
\seealso{
	 \code{\link{mesoBin}}.
}
