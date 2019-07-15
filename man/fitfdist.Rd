\name{fitFDist}
\alias{fitFDist}
\alias{fitFDistRobustly}
\title{Moment Estimation of Scaled F-Distribution}
\description{
Moment estimation of the parameters of a scaled F-distribution given one of the degrees of freedom.
This function is called internally by \code{eBayes} and \code{squeezeVar} and is not usually called directly by a user.
}
\usage{
fitFDist(x, df1, covariate=NULL)
fitFDistRobustly(x, df1, covariate=NULL, winsor.tail.p=c(0.05,0.1), trace=FALSE)
}
\arguments{
  \item{x}{numeric vector or array of positive values representing a sample from a scaled F-distribution.}
  \item{df1}{the first degrees of freedom of the F-distribution. Can be a single value, or else a vector of the same length as \code{x}.}
  \item{covariate}{if non-\code{NULL}, the estimated scale value will depend on this numeric covariate.}
  \item{winsor.tail.p}{numeric vector of length 1 or 2, giving left and right tail proportions of \code{x} to Winsorize.}
  \item{trace}{logical value indicating whether a trace of the iteration progress should be printed.}
}
\details{
\code{fitFDist} implements an algorithm proposed by Smyth (2004) and Phipson et al (2016).
It estimates \code{scale} and \code{df2} under the assumption that \code{x} is distributed as \code{scale} times an F-distributed random variable on \code{df1} and \code{df2} degrees of freedom.
The parameters are estimated using the method of moments, specifically from the mean and variance of the \code{x} values on the log-scale.

When \code{covariate} is supplied, a spline curve trend will be estimated for the \code{x} values and the estimation will be adjusted for this trend (Phipson et al, 2016).

\code{fitFDistRobustly} is similar to \code{fitFDist} except that it computes the moments of the Winsorized values of \code{x}, making it robust against left and right outliers.
Larger values for \code{winsor.tail.p} produce more robustness but less efficiency.
When \code{covariate} is supplied, a loess trend is estimated for the \code{x} values.
The robust method is described by Phipson et al (2016).

As well as estimating the F-distribution for the bulk of the cases, i.e., with outliers discounted, \code{fitFDistRobustly} also returns an estimated F-distribution with reduced df2 that might be appropriate for each outlier case.
}
\note{
The algorithm used by \code{fitFDistRobustly} was revised slightly in limma 3.27.6.
The \code{prob.outlier} value, which is the lower bound for \code{df2.shrunk}, may be slightly smaller than previously.
}
\value{
\code{fitFDist} produces a list with the following components:
  \item{scale}{scale factor for F-distribution. A vector if \code{covariate} is non-\code{NULL}, otherwise a scalar.}
  \item{df2}{the second degrees of freedom of the fitted F-distribution.}

\code{fitFDistRobustly} returns the following components as well:
  \item{tail.p.value}{right tail probability of the scaled F-distribution for each \code{x} value.}
  \item{prob.outlier}{posterior probability that each case is an outlier relative to the scaled F-distribution with degrees of freedom \code{df1} and \code{df2}.}
  \item{df2.outlier}{the second degrees of freedom associated with extreme outlier cases.}
  \item{df2.shrunk}{numeric vector of values for the second degrees of freedom, with shrunk values for outliers. Most values are equal to \code{df2}, but outliers have reduced values depending on how extreme each case is. All values lie between \code{df2.outlier} and \code{df2}.}
}

\author{Gordon Smyth and Belinda Phipson}

\references{
Smyth, G. K. (2004). Linear models and empirical Bayes methods for
     assessing differential expression in microarray experiments.
     \emph{Statistical Applications in Genetics and Molecular Biology}, \bold{3},
     No. 1, Article 3.
     \url{http://www.statsci.org/smyth/pubs/ebayes.pdf}

Phipson, B, Lee, S, Majewski, IJ, Alexander, WS, and Smyth, GK (2016).
Robust hyperparameter estimation protects against hypervariable genes and improves power to detect differential expression.
\emph{Annals of Applied Statistics} 10, 946-963.
\url{http://projecteuclid.org/euclid.aoas/1469199900}
}

\seealso{
This function is called by \code{\link{squeezeVar}}, which in turn is called by \code{\link{eBayes}} and \code{\link{treat}}.

This function calls \code{\link{trigammaInverse}}.
}
\examples{
x <- rf(100,df1=8,df2=16)
fitFDist(x,df1=8)
}
\keyword{distribution}
