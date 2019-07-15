\name{TestResults-class}
\docType{class}
\alias{TestResults-class}
\alias{show,TestResults-method}
\alias{summary.TestResults}
\alias{levels.TestResults}
\alias{labels.TestResults}
\title{Matrix of Test Results - class}

\description{
A matrix-based class for storing the results of simultanous tests.
\code{TestResults} objects are usually created by \code{\link{decideTests}}.
}

\usage{
\S3method{summary}{TestResults}(object, \dots)
\S3method{labels}{TestResults}(object, \dots)
\S3method{levels}{TestResults}(x)
}

\arguments{
\item{object, x}{object of class \code{TestResults}}
\item{\dots}{other arguments are not used}
}

\section{Slots/List Components}{
A \code{TestResults} object is essentially a numeric matrix with elements equal to \code{0}, \code{1} or \code{-1}.
Zero represents acceptance of the null hypothesis, \code{1} indicates rejection in favor of the right tail alternative and \code{-1} indicates rejection in favor of the left tail alternative.

\code{TestResults} objects can be created by \code{new("TestResults",results)} where \code{results} is a matrix.
Objects of this class contain no slots (other than \code{.Data}), although the attributes \code{dim} and \code{dimnames} may be treated as slots.
}

\section{Methods}{
This class inherits directly from class \code{matrix} so any operation appropriate for matrices will work on objects of this class.
\code{\link[methods]{show}} and \code{summary} methods are also implemented.

Functions in LIMMA which operate on \code{TestResults} objects include
\code{\link{heatDiagram}},
\code{\link{vennCounts}},
\code{\link{vennDiagram}},
\code{\link{write.fit}}.
}

\author{Gordon Smyth}

\seealso{
  \link{02.Classes} gives an overview of all the classes defined by this package.
  \link{08.Tests} gives an overview of multiple testing.
}

\examples{\dontrun{
#  Assume a data object y and a design matrix
fit <- lmFit(y, design)
fit <- eBayes(fit)
results <- decideTests(fit)
summary(results)
}}

\keyword{classes}
\keyword{htest}
