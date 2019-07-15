\name{arrayWeightsQuick}
\alias{arrayWeightsQuick}
\title{Array Quality Weights}
\description{
Estimates relative quality weights for each array in a multi-array experiment with replication.
}
\usage{
arrayWeightsQuick(y, fit)
}
\arguments{
  \item{y}{the data object used to estimate \code{fit}.
  Can be of any class which can be coerced to matrix, including \code{matrix}, \code{MAList}, \code{marrayNorm} or \code{ExpressionSet}.}
  \item{fit}{\code{MArrayLM} fitted model object}
}
\details{
Estimates the relative reliability of each array by measuring how well the expression values for that array follow the linear model.

This is a quick and dirty version of \code{\link{arrayWeights}}.
}
\value{
Numeric vector of weights of length \code{ncol(fit)}.
}
\references{
Ritchie, M. E., Diyagama, D., Neilson, van Laar, R., J., Dobrovic, A., Holloway, A., and Smyth, G. K. (2006).
Empirical array quality weights in the analysis of microarray data.
BMC Bioinformatics 7, 261.
\url{http://www.biomedcentral.com/1471-2105/7/261}
}
\author{Gordon Smyth}
\seealso{
See \link{arrayWeights}.
An overview of LIMMA functions for reading data is given in \link{03.ReadingData}.
}
\examples{
\dontrun{
fit <- lmFit(y, design)
arrayWeightsQuick(y, fit)
}
}
