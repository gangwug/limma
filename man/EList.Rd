\name{EList-class}
\docType{class}
\alias{EList-class}
\alias{EListRaw-class}
\title{Expression List - class}

\description{
A list-based S4 classes for storing expression values (E-values), for example for a set of one-channel microarrays or a set of RNA-seq samples.
\code{EListRaw} holds expression values on the raw scale.
\code{EList} holds expression values on the log scale, usually after background correction and normalization.

\code{EListRaw} objects are often created by \code{\link{read.maimages}}, while
\code{EList} objects are often created by \code{\link{normalizeBetweenArrays}} or by \code{\link{voom}}.
Alternatively, an \code{EList} object can be created directly by \code{new("EList",x)}, where \code{x} is a list.
}

\section{Required Components}{
These classes contains no slots (other than \code{.Data}), but objects should contain a list component \code{E}:
\describe{
  \item{\code{E}}{numeric matrix containing expression values.
  In an \code{EListRaw} object, the expression values are unlogged, while in an \code{EList} object, they are log2 values.
  Rows correspond to probes and columns to samples.}
}
}

\section{Optional Components}{
Optional components include:
\describe{
  \item{\code{Eb}}{numeric matrix containing unlogged background expression values, of same dimensions as \code{E}. For an \code{EListRaw} object only.}
  \item{\code{weights}}{numeric matrix of same dimensions as \code{E} containing relative spot quality weights.  Elements should be non-negative.}
  \item{\code{other}}{list containing other matrices, all of the same dimensions as \code{E}.}
  \item{\code{genes}}{data.frame containing probe information. Should have one row for each probe. May have any number of columns.}
  \item{\code{targets}}{data.frame containing information on the target RNA samples.  Rows correspond to samples.  May have any number of columns.}
}

Valid \code{EList} or \code{EListRaw} objects may contain other optional components, but all probe or sample information should be contained in the above components.
}

\section{Methods}{
These classes inherit directly from class \code{list} so any operation appropriate for lists will work on objects of this class.
In addition, \code{EList} objects can be \link[limma:subsetting]{subsetted} and \link[limma:cbind]{combined}.
\code{EList} objects will return dimensions and hence functions such as \code{\link[limma:dim]{dim}}, \code{\link[base:nrow]{nrow}} and \code{\link[base:nrow]{ncol}} are defined. 
\code{EList}s also inherit a \code{\link[methods]{show}} method from the virtual class \code{\link[limma:LargeDataObject]{LargeDataObject}}, which means that \code{ELists} will print in a compact way.
}

\author{Gordon Smyth}

\seealso{
  \link{02.Classes} gives an overview of all the classes defined by this package.
  
  \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}} is a more formal class in the Biobase package used for the same purpose.
}

\keyword{classes}

