\name{logsumexp}
\alias{logsumexp}
\title{Log Sum of Exponentials}
\description{Compute \code{log( exp(x)+exp(y) )} without floating overflow or underflow}
\usage{
logsumexp(x, y)
}
\arguments{
  \item{x}{a numeric vector or matrix.}
  \item{y}{a numeric vector or matrix of same size as \code{x}.}
}

\details{
The computation uses \code{logcosh()}.
}

\value{Numeric vector or matrix of same dimensions as \code{x}.}

\author{Gordon K Smyth}

\seealso{
\code{\link{logcosh}}
}

\examples{
x <- y <- c(1e-8,1e-7,1e-6,1e-5,1e-4,1,3,50,800)
logsumexp(x,y)
log( exp(x)+exp(y) )
}

