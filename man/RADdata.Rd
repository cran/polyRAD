\name{RADdata}
\alias{RADdata}
\alias{plot.RADdata}
\title{
RADdata object constructor
}
\description{
\code{RADdata} is used internally to generate objects of the S3 class 
\dQuote{RADdata} by \pkg{polyRAD} functions for importing read depth data.
It is also available at the user level for cases where the data for import are
not already in a format supported by \pkg{polyRAD}.
}
\usage{
RADdata(alleleDepth, alleles2loc, locTable, possiblePloidies, contamRate,
        alleleNucleotides, taxaPloidy)
        
\method{plot}{RADdata}(x, \dots)
}
\arguments{
  \item{alleleDepth}{
An integer matrix, with taxa in rows and alleles in columns.  Taxa names should
be included as row names.  Each value indicates the number of reads for a given
allele in a given taxon.  There should be no \code{NA} values; use zero to 
indicate no reads.
}
  \item{alleles2loc}{
An integer vector with one value for each column of \code{alleleDepth}.  The 
number indicates the identity of the locus to which the allele belongs.  A 
locus can have any number of alleles assigned to it (including zero).
}
  \item{locTable}{
A data frame, where locus names are row names.  There must be at least as 
many rows as the highest value of \code{alleles2loc}; each number in 
\code{alleles2loc} corresponds to a row index in \code{locTable}.  No columns
are required, although if provided a column named \dQuote{Chr} will be used for
indicating chromosome identities, a column named \dQuote{Pos} will be used
for indicating physical position, and a column named \dQuote{Ref} will be used
to indicate the reference sequence.
}
  \item{possiblePloidies}{
A list, where each item in the list is an integer vector (or a numeric vector
that can be converted to integer).  Each vector indicates an inheritance 
pattern that markers in the dataset might obey.  \code{2} indicates diploid,
\code{4} indicates autotetraploid, \code{c(2, 2)} indicates allotetraploid,
\emph{etc.}
}
  \item{contamRate}{
A number ranging from zero to one (although in practice probably less than 
0.01) indicating the expected sample cross-contamination rate.
}
  \item{alleleNucleotides}{
A character vector with one value for each column of 
\code{alleleDepth}, indicating the DNA sequence for that allele.  Typically 
only the sequence at variable sites is provided, although intervening 
non-variable sequence can also be provided.
}
  \item{taxaPloidy}{
An integer vector indicating ploidies of taxa.  If a single value is provided,
it will be assumed that all taxa are the same ploidy.  Otherwise, one value
must be provided for each taxon.  If unnamed, it is assumed that taxa are in
the same order as the rows of \code{alleleDepth}.  If named, names must match
the row names of \code{alleleDepth} but do not need to be in the same order.
This value is used as a multiplier with \code{possiblePloidies}; see Details.
}
\item{x}{
A \dQuote{RADdata} object.
}
\item{\dots}{
Additional arguments to pass to \code{plot}, for example \code{col} or 
\code{pch}.
}
}

\details{
For a single locus, ideally the string provided in \code{locTable$Ref} and all
strings in \code{alleleNucleotides} are the same length, so that SNPs and indels
may be matched by position.  The character \dQuote{-} indicates a deletion with
respect to the reference, and can be used within \code{alleleNucleotides}.  The
character \dQuote{.} is a placeholder where other alleles have an insertion with
respect to the reference, and may be used in \code{locTable$Ref} and
\code{alleleNucleotides}.  Note that it is possible for the sequence in
\code{locTable$Ref} to be absent from \code{alleleNucleotides} if the reference
haplotype is absent from the dataset, as may occur if the reference genome is that
of a related species and not the actual study species.  For the
\code{alleleNucleotides} vector, the attribute \code{"Variable_sites_only"}
indicates whether non-variable sequence in between variants is included; this
needs to be \code{FALSE} for other functions to determine the position of each
variant within the set of tags.

Inheritance mode is determined by multiplying the values in
\code{possiblePloidies} by the values in \code{taxaPloidy} and dividing by two.
For example, if you wanted to assume autotetraploid inheritance across the
entire dataset, you could set \code{possiblePloidies = list(4)} and
\code{taxaPloidy = 2}, or alternatively \code{possiblePloidies = list(2)} and
\code{taxaPloidy = 4}.  To indicate a mix of diploid and allotetraploid
inheritance across loci, set \code{possiblePloidies = list(2, c(2, 2))} and
\code{taxaPloidy = 2}.  If taxa themselves vary in ploidy, provide one
value of \code{taxaPloidy} for each taxon.  All inheritance modes listed in
\code{possiblePloidies} apply equally to all taxa, even when ploidy varies
by taxon.
}

\value{
An object of the S3 class \dQuote{RADdata}.  The following slots are available
using the \code{$} operator:

\item{alleleDepth }{Identical to the argument provided to the function.}
\item{alleles2loc }{Identical to the argument provided to the function.}
\item{locTable }{Identical to the argument provided to the function.}
\item{possiblePloidies }{The \code{possiblePloidies} argument, converted to
integer.}
\item{locDepth }{A matrix with taxa in rows and loci in columns, with read 
depth summed across all alleles for each locus.  Column names are locus 
numbers rather than locus names.  See \code{\link{GetLocDepth}} for retrieving
the same matrix but with locus names as column names.}
\item{depthSamplingPermutations }{A numeric matrix with taxa in rows and 
alleles in columns.  It is calculated as \eqn{log(locDepth choose alleleDepth)}.
This is used as a coefficient for likelihood estimations done by other
\pkg{polyRAD} functions (\emph{i.e.} \code{\link{AddGenotypeLikelihood}}).}
\item{depthRatio }{A numeric matrix with taxa in rows and alleles in columns.
Calculated as \eqn{alleleDepth / locDepth}.  Used by other \pkg{polyRAD}
functions for rough estimation of genotypes and allele frequency.}
\item{antiAlleleDepth }{An integer matrix with taxa in rows and alleles in
columns.  For each allele, the number of reads from the locus that do NOT
belong to that allele.  Calculated as \eqn{locDepth - alleleDepth}. Used for
likelihood estimations by other \pkg{polyRAD} functions.}
\item{alleleNucleotides }{Identical to the argument provided to the 
function.}
\item{taxaPloidy }{A named integer vector with one value per taxon, indicating
the ploidy of taxa.}

The object additionally has several attributes (see \code{\link{attr}}):

\item{taxa }{A character vector listing all taxa names, in the same order as
the rows of \code{alleleDepth}.}
\item{nTaxa }{An integer indicating the number of taxa.}
\item{nLoc }{An integer indicating the number of loci in \code{locTable}.}
\item{contamRate }{Identical to the argument provided to the function.}

The \code{plot} method performs a principal components analysis with 
\code{\link{AddPCA}} if not already done, then plots the first two axes.
Points represent individuals (taxa).  If mapping population parents have been
noted in the object (see \code{\link{SetDonorParent}}), they are indicated in
the plot.
}
\author{
Lindsay V. Clark
}

\seealso{
Data import functions that internally call \code{RADdata}:

\code{\link{readHMC}}, \code{\link{readTagDigger}}, 
\code{\link{VCF2RADdata}}, \code{\link{readStacks}},
\code{\link{readTASSELGBSv2}}, \code{\link{readProcessSamMulti}},
\code{\link{readProcessIsoloci}}
}
\examples{
# create the dataset
mydepth <- matrix(sample(100, 16), nrow = 4, ncol = 4,
                  dimnames = list(paste("taxon", 1:4, sep = ""),
                  paste("loc", c(1,1,2,2), "_", c(0,1,0,1), sep = "")))
mydata <- RADdata(mydepth, c(1L,1L,2L,2L), 
                  data.frame(row.names = c("loc1", "loc2"), Chr = c(1,1),
                             Pos = c(2000456, 5479880)),
                  list(2, c(2,2)), 0.001, c("A", "G", "G", "T"), 6)

# inspect the dataset
mydata
mydata$alleleDepth
mydata$locDepth
mydata$depthRatio
mydata$taxaPloidy

# the S3 class structure is flexible; other data can be added
mydata$GPS <- data.frame(row.names = attr(mydata, "taxa"),
                         Lat = c(43.12, 43.40, 43.05, 43.27),
                         Long = -c(70.85, 70.77, 70.91, 70.95))
mydata$GPS

# If you have NA in your alleleDepth matrix to indicate zero reads,
# perform the following before running the RADdata constructor:
mydepth[is.na(mydepth)] <- 0L

# plotting a RADdata object
plot(mydata)
}

\keyword{ methods }

