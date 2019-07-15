\name{avereps}
\alias{avereps}
\alias{avereps.default}
\alias{avereps.MAList}
\alias{avereps.EList}
\alias{avereps.EListRaw}
\alias{avereps.RGList}
\title{Average Over Irregular Replicate Probes}
\description{
Condense a microarray data object so that values for within-array replicate probes are replaced with their average.
}
\usage{
\method{avereps}{default}(x, ID=rownames(x), \dots)
\method{avereps}{MAList}(x, ID=NULL, \dots)
\method{avereps}{EList}(x, ID=NULL, \dots)
}
\arguments{
  \item{x}{a matrix-like object, usually a matrix, \code{MAList} or \code{EList} object.}
  \item{ID}{probe identifier.}
  \item{\dots}{other arguments are not currently used.}
}
\details{
A new data object is computed in which each probe ID is represented by the average of its replicate spots or features.

For an \code{MAList} object, the components \code{M} and \code{A} are both averaged in this way, as are \code{weights} and any matrices found in \code{object$other}.
For an \code{MAList} object, \code{ID} defaults to \code{MA$genes$ID} is that exists, otherwise to \code{rownames(MA$M)}.

\code{EList} objects are similar, except that the \code{E} component is averaged instead of \code{M} and \code{A}.

If \code{x} is of mode \code{"character"}, then the replicate values are assumed to be equal and the first is taken as the average.
}

\note{
This function should only be applied to normalized log-expression values, and not to raw unlogged expression values.
It will generate an error message if applied to \code{RGList} or \code{EListRaw} objects.
}

\value{
A data object of the same class as \code{x} with a row for each unique value of \code{ID}.
}

\author{Gordon Smyth}
\seealso{
  \code{\link{avedups}}, \code{\link{avearrays}}. Also \code{\link[base]{rowsum}} in the base package.
  
  \link{02.Classes} gives an overview of data classes used in LIMMA.
}

\examples{
x <- matrix(rnorm(8*3),8,3)
colnames(x) <- c("S1","S2","S3")
rownames(x) <- c("b","a","a","c","c","b","b","b")
avereps(x)
}
\keyword{array}
