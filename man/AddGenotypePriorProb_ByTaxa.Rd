\name{AddGenotypePriorProb_ByTaxa}
\alias{AddGenotypePriorProb_ByTaxa}
\alias{AddGenotypePriorProb_ByTaxa.RADdata}
\title{
Estimate Prior Genotype Probabilities on a Per-Taxon Basis
}
\description{
Using local allele frequencies estimated by \code{\link{AddAlleleFreqByTaxa}}
and assuming Hardy-Weinberg Equilibruim or inbreeding on a local scale,
\code{AddGenotypePriorProb_ByTaxa} estimates prior genotype probabilities at
each taxon, allele, and possible ploidy.  These are then stored in the
\code{$priorProb} slot of the \code{"\link{RADdata}"} object.
}
\usage{
AddGenotypePriorProb_ByTaxa(object, ...)
\method{AddGenotypePriorProb_ByTaxa}{RADdata}(object, selfing.rate = 0, \dots)
}
\arguments{
  \item{object}{
A \code{"\link{RADdata}"} object.  \code{\link{AddAlleleFreqByTaxa}} should have
already been run.
}
\item{selfing.rate}{
A number ranging from zero to one indicating the frequency of
self-fertilization in the species.
}
  \item{\dots}{
Additional arguments (none implemented).
}
}

\value{
A \code{"RADdata"} object identical to that passed to the function, but with a list
added to the \code{$priorProb} slot.  Each item in the list corresponds to one
ploidy in \code{object$possiblePloidies}, and is a three-dimensional array with 
allele copy number in the first dimension, taxa in the second dimension, and
alleles in the third dimension.  The values in the array are prior genotype
probabilities.  Additionally, \code{object$possiblePloidies} is copied to
\code{object$priorProbPloidies}, and \code{"taxon"} is recorded in the 
\code{"priorType"} attribute.
}

\author{
Lindsay V. Clark
}

\seealso{
\code{\link{AddGenotypePriorProb_HWE}} for equations used for genotype prior
probability estimation.

\code{\link{AddGenotypePriorProb_Mapping2Parents}},
\code{\link{AddGenotypeLikelihood}}
}
\examples{
# load data
data(exampleRAD)
# do PCA
exampleRAD <- AddPCA(exampleRAD, nPcsInit = 3)
# get allele frequencies
exampleRAD <- AddAlleleFreqByTaxa(exampleRAD)

# add prior probabilities
exampleRAD <- AddGenotypePriorProb_ByTaxa(exampleRAD)

exampleRAD$priorProb[[1]][,1,]
exampleRAD$priorProb[[2]][,1,]
exampleRAD$priorProb[[1]][,2,]
exampleRAD$priorProb[[2]][,2,]

# try it with inbreeding
exampleRAD <- AddGenotypePriorProb_ByTaxa(exampleRAD, selfing.rate = 0.5)

exampleRAD$priorProb[[1]][,1,]
}

\keyword{ distribution }
\keyword{ methods }
