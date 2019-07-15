\title{Plot Expression Densities}
\name{plotDensities}
\alias{plotDensities}
\alias{plotDensities.default}
\alias{plotDensities.RGList}
\alias{plotDensities.MAList}
\alias{plotDensities.EListRaw}
\alias{plotDensities.EList}
\description{
Plot the density of expression values for multiple arrays on the same plot.
}
\usage{
\method{plotDensities}{RGList}(object, log=TRUE, group=NULL, col=NULL, main="RG Densities",
              bc.method="subtract", \dots)
\method{plotDensities}{MAList}(object, log=TRUE, group=NULL, col=NULL, main="RG Densities", \dots)
\method{plotDensities}{EListRaw}(object, log=TRUE, bc.method="subtract", \dots)
\method{plotDensities}{EList}(object, log=TRUE, \dots)
\method{plotDensities}{default}(object, group=NULL, col=NULL, main=NULL, legend="topleft", \dots)
}
\arguments{
  \item{object}{an \code{RGList}, \code{MAList}, \code{EListRaw} or \code{EList} object containing expression data.  Or any data object that can be coerced to a matrix.}
  \item{log}{logical, should densities be plotted on the log2 scale?}
  \item{group}{optional vector or factor classifying the arrays into groups.  Should be same length as \code{ncol(object)}.}
  \item{col}{optional vector of colors of the same length as the number of groups.}
  \item{main}{the main title for the plot.}
  \item{bc.method}{background subtraction method passed to \code{\link{backgroundCorrect}}.}
  \item{legend}{character string giving position to place legend.
  See \code{\link{legend}} for possible values.
  Can also be logical, with \code{FALSE} meaning no legend.}
  \item{\ldots}{other arguments are passed to \code{plotDensities.default} or \code{\link{density}}.}
}

\details{
This function is useful to display and contrast the distribution of expression values on different arrays.
It can for example be used to display the effects of between-array normalization.
See the section on between-array normalization in the LIMMA User's Guide.
}
\value{A plot is created on the current graphics device.}
\author{Natalie Thorne and Gordon Smyth}
\seealso{
An overview of diagnostic plots in LIMMA is given in \link{09.Diagnostics}.
There is a section using \code{plotDensities} in conjunction with between-array normalization
in the \link[=limmaUsersGuide]{LIMMA User's Guide}.

This function uses \code{\link{density}} and \code{\link{matplot}}.
}
\examples{
\dontrun{
#  Default is to plot red channels in red and green channels in green
plotDensities(MA)

#  Alternatively colors
plotDensities(MA, col=c("red","blue"))

#  Color by group, with three groups:
plotDensities(MA,group=group,col=c("blue","orange","green"))
}
}
\keyword{hplot}
