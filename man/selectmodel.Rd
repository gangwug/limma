\name{selectModel}
\alias{selectModel}
\title{Select Appropriate Linear Model}
\description{Select the best fitting linear model for each gene by minimizing an information criterion.}
\usage{
selectModel(y, designlist, criterion="aic", df.prior=0, s2.prior=NULL, s2.true=NULL, \dots)
}
\arguments{
  \item{y}{a matrix-like data object, containing log-ratios or log-values of expression for a series of microarrays.
  Any object class which can be coerced to matrix is acceptable including \code{numeric}, \code{matrix}, \code{MAList}, \code{marrayNorm}, \code{ExpressionSet} or \code{PLMset}.}
  \item{designlist}{list of design matrices} 
  \item{criterion}{information criterion to be used for model selection, \code{"aic"}, \code{"bic"} or \code{"mallowscp"}.}
  \item{df.prior}{prior degrees of freedom for residual variances. See \code{\link{squeezeVar}}}
  \item{s2.prior}{prior value for residual variances, to be used if \code{df.prior}>0.}
  \item{s2.true}{numeric vector of true variances, to be used if \code{criterion="mallowscp"}.}
  \item{\dots}{other optional arguments to be passed to \code{lmFit}}
}

\value{
List with components
\item{IC}{matrix of information criterion scores, rows for probes and columns for models}
\item{pref}{factor indicating the model with best (lowest) information criterion score}
}

\details{
This function chooses, for each probe, the best fitting model out of a set of alternative models represented by a list of design matrices.
Selection is by Akaike's Information Criterion (AIC), Bayesian Information Criterion (BIC) or by Mallow's Cp.

The criteria have been generalized slightly to accommodate an information prior on the variances represented by \code{s2.prior} and \code{df.prior} or by \code{s2.post}.
Suitable values for these parameters can be estimated using \code{\link{squeezeVar}}.
}

\seealso{
An overview of linear model functions in limma is given by \link{06.LinearModels}.
}

\author{Alicia Oshlack and Gordon Smyth}
\examples{
nprobes <- 100
narrays <- 5
y <- matrix(rnorm(nprobes*narrays),nprobes,narrays)
A <- c(0,0,1,1,1)
B <- c(0,1,0,1,1)
designlist <- list(
  None=cbind(Int=c(1,1,1,1,1)),
  A=cbind(Int=1,A=A),
  B=cbind(Int=1,B=B),
  Both=cbind(Int=1,AB=A*B),
  Add=cbind(Int=1,A=A,B=B),
  Full=cbind(Int=1,A=A,B=B,AB=A*B)
)
out <- selectModel(y,designlist)
table(out$pref)
}
\keyword{models}
\keyword{regression}
