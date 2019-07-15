\name{normalizeVSN}
\alias{normalizeVSN}
\alias{normalizeVSN.RGList}
\alias{normalizeVSN.EListRaw}
\alias{normalizeVSN.default}
\title{Variance Stabilizing Normalization (vsn)}

\description{
Apply variance stabilizing normalization (vsn) to limma data objects.
}

\usage{
normalizeVSN(x, \dots)
}

\arguments{
  \item{x}{a numeric \code{matrix}, \code{EListRaw} or \code{\link[limma:rglist]{RGList}} object.}
  \item{\dots}{other arguments are passed to \code{vsn}}
}

\details{
This is an interface to the \code{vsnMatrix} function from the vsn package.
The input \code{x} should contain raw intensities.
If \code{x} contains background and well as foreground intensities, these will be subtracted from the foreground intensities before \code{vsnMatrix} is called.

Note that the vsn algorithm performs background correction and normalization simultaneously.
If the data are from two-color microarrays, then the red and green intensities are treated as if they were single channel data, i.e., red and green channels from the same array are treated as unpaired.
This algorithm is therefore separate from the backgroundCorrection, normalizeWithinArrays, then normalizeBetweenArrays paradigm used elsewhere in the limma package.
}

\value{
The class of the output depends on the input.
If \code{x} is a matrix, then the result is a matrix of the same size.
If \code{x} is an \code{EListRaw} object, then an \code{EList} object with expression values on the log2 scale is produced.
For \code{x} is an \code{RGList}, then an \code{\link[limma:malist]{MAList}} object with M and A-values on the log2 scale is produced.
}

\author{Gordon Smyth}

\references{
Huber, W, von Heydebreck, A, Sueltmann, H, Poustka, A, Vingron, M (2002).
Variance stabilization applied to microarray data calibration and to the quantification of differential expression.
Bioinformatics 18 Supplement 1, S96-S104.
}

\seealso{
  An overview of LIMMA functions for normalization is given in \link{05.Normalization}.

  See also \code{\link[vsn:vsn2]{vsnMatrix}} in the vsn package.
}

\examples{
ngenes <- 100
narrays <- 4
x <- matrix(rnorm(ngenes*narrays),100,4)
y <- normalizeVSN(x)
}

\keyword{models}
\keyword{multivariate}
