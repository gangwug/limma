\name{logcosh}
\alias{logcosh}
\title{Logarithm of cosh}
\description{Compute \code{log(cosh(x))} without floating overflow or underflow}
\usage{
logcosh(x)
}
\arguments{
  \item{x}{a numeric vector or matrix.}
}

\details{
The computation uses asymptotic expressions for very large or very small arguments.
For intermediate arguments, \code{log(cosh(x))} is returned.
}

\value{Numeric vector or matrix of same dimensions as \code{x}.}

\author{Gordon K Smyth}

\seealso{
\code{\link{logsumexp}}
}

\examples{
x <- c(1e-8,1e-7,1e-6,1e-5,1e-4,1,3,50,800)
logcosh(x)
log(cosh(x))
}
