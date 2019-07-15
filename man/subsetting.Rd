\name{subsetting}
\alias{subsetting}
\alias{[.RGList}
\alias{[.MAList}
\alias{[.EListRaw}
\alias{[.EList}
\alias{[.MArrayLM}
\alias{subsetListOfArrays}
\title{Subset RGList, MAList, EListRaw, EList or MArrayLM Objects}
\description{
Return an \code{RGList}, \code{MAList}, \code{EListRaw}, \code{EList} or \code{MArrayLM} object with only selected rows and columns of the original object.
}
\usage{
\method{[}{RGList}(object, i, j)
subsetListOfArrays(object, i, j, IJ, IX, I, JX)
}
\arguments{
  \item{object}{object of class \code{RGList}, \code{MAList}, \code{EListRaw}, \code{EList} or \code{MArrayLM}.}
  \item{i,j}{elements to extract. \code{i} subsets the probes or spots while \code{j} subsets the arrays.}
  \item{IJ}{character vector giving names of components that should be subsetted by \code{i} and \code{j}.}
  \item{IX}{character vector giving names of 2-dimensional components that should be subsetted by \code{i} only.}
  \item{I}{character vector giving names of vector components that should be subsetted by \code{i}.}
  \item{JX}{character vector giving names of 2-dimensional components whose row dimension corresponds to \code{j}.}
}
\details{
\code{i,j} may take any values acceptable for the matrix components of \code{object}.
Either or both can be missing.
See the \link[base]{Extract} help entry for more details on subsetting matrices.

\code{object[]} will return the whole object unchanged.
A single index \code{object[i]} will be taken to subset rows, so \code{object[i]} and \code{object[i,]} are equivalent.

\code{subsetListOfArrays} is used internally as a utility function by the subsetting operations.
It is not intended to be called directly by users.
Values must be supplied for all arguments other than \code{i} and \code{j}.
}
\value{
An object the same as \code{object} but containing data from the specified subset of rows and columns only.
}
\author{Gordon Smyth}
\seealso{
  \code{\link[base]{Extract}} in the base package.
  
  \link{02.Classes} for a summary of the different data classes.
}
\examples{
M <- A <- matrix(11:14,4,2)
rownames(M) <- rownames(A) <- c("a","b","c","d")
colnames(M) <- colnames(A) <- c("A","B")
MA <- new("MAList",list(M=M,A=A))
MA[1:2,]
MA[c("a","b"),]
MA[1:2,2]
MA[,2]
}
