\name{mdplot}
\alias{mdplot}
\title{Mean-Difference Plot}
\description{Creates a mean-difference plot of two columns of a matrix.}

\usage{mdplot(x, columns=c(1,2), xlab="Mean", ylab="Difference", main=NULL, \dots)}

\arguments{
  \item{x}{numeric \code{matrix} with at least two columns.}
  \item{columns}{which columns of \code{x} to compare. Plot will display second minus first.}
  \item{xlab}{label for the x-axis.}
  \item{ylab}{label for the y-axis.}
  \item{main}{title of the plot. Defaults to }
  \item{\dots}{any other arguments are passed to \code{\link{plotWithHighlights}}.}
}

\details{
Plots differences vs means for a set of bivariate values.
This is a generally useful approach for comparing two correlated measures of the same underlying phenomenon.
Bland and Altman (1986) argue it is more information than a simple scatterplot of the two variables.
The bivariate values are stored as columns of \code{x}.
}

\value{A plot is created on the current graphics device.}

\references{
Cleveland, W. S., (1993). Visualizing Data. Hobart Press.

Bland, J. M., and Altman, D. G. (1986). Statistical methods for assessing agreement between two methods of clinical measurement. Lancet 327, 307-310.

See also \url{http://www.statsci.org/micrarra/refs/maplots.html}
}

\author{Gordon Smyth}

\seealso{
\code{\link{plotWithHighlights}}

\code{\link{plotMD}} is an object-oriented implementation of mean-difference plots for expression data.

An overview of diagnostic functions available in LIMMA is given in \link{09.Diagnostics}.
}

\examples{
x1 <- runif(100)
x2 <- (x1 + rnorm(100,sd=0.01))^1.2
oldpar <- par(mfrow=c(1,2))
plot(x1,x2)
mdplot(cbind(x1,x2),bg.pch=1,bg.cex=1)
par(oldpar)
}

\keyword{hplot}
