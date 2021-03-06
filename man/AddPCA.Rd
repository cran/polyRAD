\name{AddPCA}
\alias{AddPCA}
\alias{AddPCA.RADdata}
\title{
Perform Principal Components Analysis on \dQuote{RADdata} Object
}
\description{
This function uses read depth ratios or posterior genotype probabilities
(the latter preferentially) as input data for principal components analysis.
The PCA scores are then stored in the \code{$PCA} slot of the
\code{"\link{RADdata}"} object.
}
\usage{
AddPCA(object, ...)
\method{AddPCA}{RADdata}(object, nPcsInit = 10, maxR2changeratio = 0.05, 
                         minPcsOut = 1, \dots)
}
\arguments{
  \item{object}{
A \code{"\link{RADdata}"} object.
}
\item{nPcsInit}{
The number of principal component axes to initially calculate.
}
\item{maxR2changeratio}{
This number determines how many principal component axes are retained.  The 
difference in \eqn{R^2}{R-squared} values between the first and second axes 
is multiplied by \code{maxR2changeratio}.  The last axis retained is the first
axis after which the \eqn{R^2}{R-squared} value changes by less than this value.
Lower values of \code{maxR2changeratio} will result in more axes being retained.
}
\item{minPcsOut}{
The minimum number of PC axes to output, which can override 
\code{maxR2changeratio}.
}
  \item{\dots}{
Additional arguments to be passed to the \code{pca} function from the
\pkg{pcaMethods} BioConductor package.
}
}
\details{
The PPCA (probabalistic PCA) method from \pkg{pcaMethods} is used,
due to the high missing data rate that is typical of genotyping-by-sequencing
datasets.
}
\value{
A \code{"RADdata"} object identical to the one passed to the function, but with
a matrix added to the \code{$PCA} slot.  This matrix contains PCA scores, with
taxa in rows, and PC axes in columns.
}

\note{
If you see the error 

\code{Error in if (rel_ch < threshold & count > 5) \{ : 
  missing value where TRUE/FALSE needed}
  
try lowering \code{nPcsInit}.
}

\author{
Lindsay V. Clark
}

\seealso{
\code{\link{AddAlleleFreqByTaxa}}
}
\examples{
# load data
data(exampleRAD)
# do PCA
exampleRAD <- AddPCA(exampleRAD, nPcsInit = 3)

plot(exampleRAD$PCA[,1], exampleRAD$PCA[,2])
}

\keyword{ methods }

