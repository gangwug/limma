\name{as.MAList}
\alias{as.MAList}
\title{Convert marrayNorm Object to an MAList Object}

\description{
Convert marrayNorm Object to an MAList Object
}

\usage{
as.MAList(object)
}

\arguments{
\item{object}{an \code{marrayNorm} object}
}

\value{
Object of class \code{\link[=MAList-class]{MAList}}
}

\details{
The \code{marrayNorm} class is defined in the \code{marray} package.
This function converts a normalized two color microarray data object created by the \code{marray} package into the corresponding limma data object.

Note that such conversion is not necessary to access the limma linear modelling functions, because \code{lmFit} will operate on a \code{marrayNorm} data object directly.
}

\author{Gordon Smyth}

\seealso{
  \link{02.Classes} gives an overview of all the classes defined by this package.

  The marrayNorm class is defined in the marray package.
}

\keyword{classes}
\keyword{data}
