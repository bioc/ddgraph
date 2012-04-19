\name{mesoBin}
\alias{mesoBin}
\docType{data}

\title{
  A list of binary \code{DDDataSet} objects.
}
\description{
  \code{mesoBin} is a list of objects of class \code{\linkS4class{DDDataSet}}. 
  It has been generated with the following code:
  
  \code{mesoBin <- toDDDataSet(readFurlongData(), prettyNames=TRUE)}
}

\usage{data(mesoBin)}

\details{
The dataset represents binary binding signal for 5 transcription factors (TFs) at 1-5 time points
during embryonic mesoderm development in Drosophila Melanogaster (Zinzen et al, 2009). The original data 
has been binarized by taking any signal greater than the threshold authors used as positive binding event. 

The list contains 7 objects of type \code{\linkS4class{DDDataSet}} for 7 cis-regulatory
module (CRM) classes. These classes are: neg (negative class of CRMs), Meso (CRMs active in
early mesoderm), Meso_SM (CRMs active in early mesoderm and somatic muscle), VM (visceral 
muscle), SM (somatic muscle), VM_SM (active in both somatic and visceral muscle) and CM
(active in cardiac muscle). 
}

\examples{
data(mesoBin)
names(mesoBin)
class(mesoBin$VM)
}
\references{Robert P. Zinzen et al., "Combinatorial binding predicts spatio-temporal cis-regulatory activity,"
Nature 462, no. 7269 (November 5, 2009): 65-70.}
\keyword{datasets}
\seealso{
	 \code{\link{mesoCont}}.
}
