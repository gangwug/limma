\name{getEAWP}
\alias{getEAWP}
\title{Extract Basic Data from Expression Data Objects}
\description{
Given an expression data object of any known class, get the expression values, weights, probe annotation and A-values that are needed for linear modelling.
This function is called by the linear modelling functions in LIMMA.
}
\usage{
getEAWP(object)
}
\arguments{
  \item{object}{any matrix-like object containing log-expression values.
  Can be an object of class \code{MAList}, \code{EList}, \code{marrayNorm}, \code{PLMset}, \code{vsn}, or any class inheriting from \code{ExpressionSet}, or any object that can be coerced to a numeric matrix.}
}
\details{
Rows correspond to probes and columns to RNA samples.

In the case of two-color microarray data objects (\code{MAList} or \code{marrayNorm}), \code{Amean} is the vector of row means of the matrix of A-values.
For other data objects, \code{Amean} is the vector of row means of the matrix of expression values.

From April 2013, the rownames of the output \code{exprs} matrix are required to be unique.
If \code{object} has no row names, then the output rownames of \code{exprs} are \code{1:nrow(object)}.
If \code{object} has row names but with duplicated names, then the rownames of \code{exprs} are set to \code{1:nrow(object)} and the original row names are preserved in the \code{ID} column of \code{probes}.

\code{object} should be a normalized data object.
\code{getEAWP} will return an error if \code{object} is a non-normalized data object such as \code{RGList} or \code{EListRaw}, because these do not contain log-expression values.
}
\value{
A list with components
\item{exprs}{numeric matrix of log-ratios, log-intensities or log-expression values}
\item{weights}{numeric matrix of weights}
\item{probes}{data.frame of probe-annotation}
\item{Amean}{numeric vector of average log-expression for each probe}
\code{exprs} is the only required component.
The other components will be \code{NULL} if not found in the input object.
}
\author{Gordon Smyth}
\seealso{
  \link{02.Classes} gives an overview of data classes used in LIMMA.
}
