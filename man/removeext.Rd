\name{removeExt}
\alias{removeExt}

\title{Remove Common Extension from File Names}
\description{Finds and removes any common extension from a vector of file names.}

\usage{
removeExt(x, sep=".")
}

\arguments{
\item{x}{character vector}
\item{sep}{character string that separates the body of each character string from the extension.}
}

\details{
This function is used for simplifying file names, or any vector of character strings, when the strings all finish with the same suffix or extension.
If the same extension is not shared by every element of \code{x}, then it is not removed from any element.

Note that \code{sep} is interpreted as a literal character string: it is not a regular expression.
}

\value{
A character vector of the same length as \code{x} in which any common extension has been stripped off.
}

\seealso{
An overview of LIMMA functions for reading data is given in \link{03.ReadingData}.
}

\examples{
x <- c("slide1.spot","slide2.spot","slide3.spot")
removeExt(x)

x <- c("Harry - a name from Harry Potter","Hermione - a name from Harry Potter")
removeExt(x, sep=" - ")
}

\author{Gordon Smyth}

\keyword{character}
\keyword{file}
