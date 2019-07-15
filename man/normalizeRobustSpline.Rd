\name{normalizeRobustSpline}
\alias{normalizeRobustSpline}
\title{Normalize Single Microarray Using Shrunk Robust Splines}

\description{
Normalize the M-values for a single microarray using robustly fitted regression splines and empirical Bayes shrinkage.
}

\usage{
normalizeRobustSpline(M,A,layout=NULL,df=5,method="M")
}

\arguments{
  \item{M}{numeric vector of M-values}
  \item{A}{numeric vector of A-values}
  \item{layout}{list specifying the dimensions of the spot matrix and the grid matrix. Defaults to a single group for the whole array.}
  \item{df}{degrees of freedom for regression spline, i.e., the number of regression coefficients and the number of knots}
  \item{method}{choices are \code{"M"} for M-estimation or \code{"MM"} for high breakdown point regression}
}

\details{
This function implements an idea similar to print-tip loess normalization but uses regression splines in place of the loess curves and uses empirical Bayes ideas to shrink the individual print-tip curves towards a common value.
This allows the technique to introduce less noise into good quality arrays with little spatial variation while still giving good results on arrays with strong spatial variation.

The original motivation for the robustspline method was to use whole-array information to moderate the normalization curves used for the individual print-tip groups.
This was an important issue for academically printed spotted two-color microarrays, especially when some of the print-tip groups contained relatively few spots.
In these situations, robust spline normalization ensures stable results even for print-tip groups with few spots.

Modern commercial two colour arrays do not usually have print tips, so in effect the whole array is a single print-tip group, and so the need for moderating individual curves is gone.
Robustspline normalization can still be used for data from these arrays, in which case a single normalization curve is estimated.
In this situation, the method is closely analogous to global loess, with a regression spline replacing the loess curve and with robust
regression replacing the loess robustifying weights.
Robust spline normalization with \code{method="MM"} has potential advantages over global loess normalization when there a lot of differential expression or the differential expression is assymetric, because of the increased level of robustness.
The potential advantages of this approach have not been fully explored in a refereed publication however.
}

\value{
Numeric vector containing normalized M-values.
}

\author{Gordon Smyth}

\references{
Ritchie, ME, Phipson, B, Wu, D, Hu, Y, Law, CW, Shi, W, and Smyth, GK (2015).
limma powers differential expression analyses for RNA-sequencing and microarray studies.
\emph{Nucleic Acids Research} 43, e47.
\url{http://nar.oxfordjournals.org/content/43/7/e47}
}

\seealso{
\code{normalizeRobustSpline} uses \code{ns} in the splines package to specify regression splines and \code{rlm} in the MASS package for robust regression.

This function is usually accessed through \code{\link{normalizeWithinArrays}}.
An overview of LIMMA functions for normalization is given in \link{05.Normalization}.
}

\examples{
A <- 1:100
M <- rnorm(100)
normalized.M <- normalizeRobustSpline(M,A)

# Usual usage
\dontrun{MA <- normalizeWithinArrays(RG, method="robustspline")}
}

\keyword{normalization}
