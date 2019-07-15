\name{fitGammaIntercept}
\alias{fitGammaIntercept}
\title{Fit Intercept to Vector of Gamma Distributed Variates}
\description{
Fit Intercept to Vector of Gamma Distributed Variates}
\usage{
fitGammaIntercept(y,offset=0,maxit=1000)
}
\arguments{
  \item{y}{numeric vector of positive response values.}
  \item{offset}{numeric vector giving known part of the expected value of \code{y}. Can be a single value, or else a vector of the same length as \code{y}.}
  \item{maxit}{maximum number of Newton iterations to be done.}
}
\details{
The values \code{y} are assumed to follow a gamma distribution with common shape parameter and with expected values given by \code{x+offset}.
The function implements a globally convergent Newton iteration to estimate \code{x}.
}
\value{
Numeric value giving intercept.
}
\author{Gordon Smyth and Belinda Phipson}

\references{
Phipson, B. (2013).
\emph{Empirical Bayes modelling of expression profiles and their associations}.
PhD Thesis. University of Melbourne, Australia.
}

\seealso{
This function is called by \code{\link{genas}}.
}
\examples{
offset <- runif(10)
x <- 9
mu <- x+offset
y <- rgamma(10,shape=20,scale=mu/20)
fitGammaIntercept(y,offset=offset)
}
\keyword{distribution}
