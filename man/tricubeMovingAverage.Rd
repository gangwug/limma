\name{tricubeMovingAverage}
\alias{tricubeMovingAverage}
\title{Moving Average Smoother With Tricube Weights}
\description{
Apply a moving average smoother with tricube distance weights to a numeric vector.
}
\usage{
tricubeMovingAverage(x, span=0.5, power=3)
}
\arguments{
  \item{x}{numeric vector}
  \item{span}{the smoother span. This gives the proportion of \code{x} values that contribute to each moving average. Larger values give more smoothness. Should be positive but not greater than 1.}
  \item{power}{a positive exponent used to compute the tricube weights. \code{power=3} gives the usual tricube weights.  Smaller values give more even weighting.  Should be greater than 0.}
}
\details{
This function smooths a vector (considered as a time series) using a moving average with tricube weights.
Specifically, the function computes running weighted means of \code{w} consecutive values of \code{x}, where the window width \code{w} is equal to \code{2*h+1} with \code{h = 2*floor(span*length(x)/2)}.
The window width \code{w} is always odd so that each window has one of the original \code{x} values at its center.
Each weighted mean uses a set of tricube weights so that values near the ends of the window receive less weight.

The smoother returns a vector of the same length as input.
At the start and end of the vector, the series is considered to be extended by missing values, and the weighted average is computed only over the observed values.
In other words, the window width is reduced to \code{h+1} at the boundaries with asymmetric weights.

The result of this function is similar to a least squares loess curve of degree zero, with a couple of differences.
First, a continuity correction is applied when computing the distance to neighbouring points, so that exactly \code{w} points are included with positive weights in each average.
Second, the span halves at the end points so that the smoother is more sensitive to trends at the ends.

The \code{filter} function in the stats package is called to do the low-level calculations.

This function is used by \code{\link{barcodeplot}} to compute enrichment worms.
}
\value{
Numeric vector of same length as \code{x} containing smoothed values.
}
\examples{
x <- rbinom(100,size=1,prob=0.5)
plot(1:100,tricubeMovingAverage(x))
}
\seealso{
\code{\link{filter}}, \code{\link{barcodeplot}}, \code{\link[edgeR:loessByCol]{loessByCol}}
}
\author{Gordon Smyth}
\keyword{smooth}
\concept{Loess curve}

