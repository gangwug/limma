\name{tmixture}
\alias{tmixture.vector}
\alias{tmixture.matrix}
\title{Estimate Scale Factor in Mixture of t-Distributions}
\description{
These functions estimate the unscaled standard deviation of the true (unobserved) log fold changes for differentially expressed genes.
They are used by the functions \code{ebayes} and \code{eBayes} and are not intended to be called directly by users.
}
\usage{
tmixture.vector(tstat, stdev.unscaled, df, proportion, v0.lim = NULL)
tmixture.matrix(tstat, stdev.unscaled, df, proportion, v0.lim = NULL)
}
\arguments{
  \item{tstat}{numeric vector or matrix of t-statistics.  \code{tmixture.vector} assumes a vector while \code{tmixture.matrix} assumes a matrix.}
  \item{stdev.unscaled}{numeric vector or matrix, conformal with \code{tstat}, containing the unscaled standard deviations of the coefficients used to compute the t-statistics.}
  \item{df}{numeric vector giving the degrees of freedom associated with \code{tstat}.}
  \item{proportion}{assumed proportion of genes that are differentially expressed.}
  \item{v0.lim}{numeric vector of length 2 giving the lower and upper limits for the estimated unscaled standard deviations.}
}
\details{
The values in each column of \code{tstat} are assumed to follow a mixture of an ordinary t-distribution, with mixing proportion \code{1-proportion}, and \code{(v0+v1)/v1} times a t-distribution, with mixing proportion \code{proportion}.
Here \code{v1} is \code{stdev.unscaled^2} and \code{v0} is the value to be estimated.
}
\value{
Numeric vector, of length equal to the number of columns of \code{tstat}, containing estimated \code{v0} values.
}
\seealso{
\code{\link{eBayes}}
}
\author{Gordon Smyth}
\keyword{distribution}
