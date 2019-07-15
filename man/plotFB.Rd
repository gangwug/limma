\title{FB-Plot}
\name{plotFB}
\alias{plotFB}
\alias{plotFB.RGList}
\alias{plotFB.EListRaw}
\alias{plotFB.default}
\description{
Creates foreground-background plots.
}
\usage{
\method{plotFB}{RGList}(x, array=1, lim="separate", pch=16, cex=0.2, \dots)
\method{plotFB}{EListRaw}(x, array=1, pch=16, cex=0.2, \dots)
}
\arguments{
  \item{x}{an \code{RGList} or \code{EListRaw} object.}
  \item{array}{integer giving the array to be plotted.}
  \item{lim}{character string indicating whether the red and green plots should have \code{"separate"} or \code{"common"} x- and y- co-ordinate limits.}
  \item{pch}{vector or list of plotting characters. Defaults to integer code 16.}
  \item{cex}{numeric vector of plot symbol expansions.} 
  \item{\dots}{any other arguments are passed to \code{plot}}
}

\details{
A foreground-background plot is a plot of log2-foreground vs log2-background for a particular array.
For two-color arrays, this function produces a pair of plots, one for the green channel and one for the red.

See \code{\link[graphics]{points}} for possible values for \code{pch}, \code{col} and \code{cex}.
}

\value{A plot is created on the current graphics device.}
\author{Gordon Smyth}
\seealso{
An overview of diagnostic functions available in LIMMA is given in \link{09.Diagnostics}.
}
\keyword{hplot}
