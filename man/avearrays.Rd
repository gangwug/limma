\name{avearrays}
\alias{avearrays}
\alias{avearrays.default}
\alias{avearrays.MAList}
\alias{avearrays.EList}
\title{Average Over Replicate Arrays}
\description{
Condense a microarray data object so that technical replicate arrays are replaced with (weighted) averages.
}
\usage{
\method{avearrays}{default}(x, ID=colnames(x), weights=NULL)
\method{avearrays}{MAList}(x, ID=colnames(x), weights=x$weights)
\method{avearrays}{EList}(x, ID=colnames(x), weights=x$weights)
}
\arguments{
  \item{x}{a matrix-like object, usually a matrix, \code{MAList} or \code{EList} object.}
  \item{ID}{sample identifier.}
  \item{weights}{numeric matrix of non-negative weights}
}
\details{
A new data object is computed in which technical replicate arrays are replaced by their (weighted) averages.

For an \code{MAList} object, the components \code{M} and \code{A} are both averaged in this way, as are \code{weights} and any matrices found in \code{object$other}.

\code{EList} objects are similar, except that the \code{E} component is averaged instead of \code{M} and \code{A}.

If \code{x} is of mode \code{"character"}, then the replicate values are assumed to be equal and the first is taken as the average.
}
\value{
A data object of the same class as \code{x} with a column for each unique value of \code{ID}.
}
\author{Gordon Smyth}
\seealso{
  \code{\link{avereps}}.
  
  \link{02.Classes} gives an overview of data classes used in LIMMA.
}

\examples{
x <- matrix(rnorm(8*3),8,3)
colnames(x) <- c("a","a","b")
avearrays(x)
}
