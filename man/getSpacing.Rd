\name{getSpacing}
\alias{getSpacing}
\title{Get Numerical Spacing}
\description{
Convert character to numerical spacing measure for within-array replicate spots.
}
\usage{
getSpacing(spacing, layout)
}
\arguments{
  \item{spacing}{character string or integer.
  Acceptable character strings are \code{"columns"}, \code{"rows"}, \code{"subarrays"} or \code{"topbottom"}.
  Integer values are simply passed through.}
  \item{layout}{list containing printer layout information}
}
\details{
\code{"rows"} means that duplicate spots are printed side-by-side by rows.
These will be recorded in consecutive rows in the data object.

\code{"columns"} means that duplicate spots are printed side-by-sidy by columns.
These will be separated in the data object by \code{layout$nspot.r} rows.

\code{"subarrays"} means that a number of sub-arrays, with identical probes in the same arrangement, are printed on each array.
The spacing therefore will be the size of a sub-array.

\code{"topbottom"} is the same as \code{"subarrays"} when there are two sub-arrays.
}
\value{
Integer giving spacing between replicate spots in the gene list.
}
\author{Gordon Smyth}
\seealso{
An overview of LIMMA functions for reading data is given in \link{03.ReadingData}.
}
\examples{
getSpacing("columns",list(ngrid.r=2,ngrid.c=2,nspot.r=20,nspot.c=19))
getSpacing("rows",list(ngrid.r=2,ngrid.c=2,nspot.r=20,nspot.c=19))
getSpacing("topbottom",list(ngrid.r=2,ngrid.c=2,nspot.r=20,nspot.c=19))
}
\keyword{IO}
