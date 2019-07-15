\name{contrasts.fit}
\alias{contrasts.fit}
\title{Compute Contrasts from Linear Model Fit}
\description{
Given a linear model fit to microarray data, compute estimated coefficients and standard errors for a given set of contrasts.
}
\usage{
contrasts.fit(fit, contrasts=NULL, coefficients=NULL)
}
\arguments{
  \item{fit}{an \code{\link[limma:marraylm]{MArrayLM}} object or a list object produced by the function \code{lm.series} or equivalent. Must contain components \code{coefficients} and \code{stdev.unscaled}.}
  \item{contrasts}{numeric matrix with rows corresponding to coefficients in \code{fit} and columns containing contrasts. May be a vector if there is only one contrast.}
  \item{coefficients}{vector indicating which coefficients are to be kept in the revised fit object. An alternative way to specify the \code{contrasts}.}
}
\value{
  An list object of the same class as \code{fit}, usually \code{\link[limma:marraylm]{MArrayLM}}. This is a list with components
  \item{coefficients}{numeric matrix containing the estimated coefficients for each contrast for each probe.}
  \item{stdev.unscaled}{numeric matrix conformal with \code{coef} containing the unscaled standard deviations for the coefficient estimators.}
  \item{\code{cov.coefficients}}{numeric \code{matrix} giving the unscaled covariance matrix of the estimable coefficients.}
Most other components found in \code{fit} are passed through unchanged, but \code{t}, \code{p.value}, \code{lods}, \code{F} and \code{F.p.value} will all be removed.
}
\details{
This function accepts input from any of the functions \code{\link{lmFit}}, \code{lm.series}, \code{mrlm}, \code{gls.series} or \code{lmscFit}.
The function re-orientates the fitted model object from the coefficients of the original design matrix to any set of contrasts of the original coefficients.
The coefficients, unscaled standard deviations and correlation matrix are re-calculated in terms of the contrasts.

The idea of this function is to fit a full-rank model using \code{lmFit} or equivalent, then use \code{contrasts.fit} to obtain coefficients and standard errors for any number of contrasts of the coefficients of the original model.
Unlike the design matrix input to \code{lmFit}, which normally has one column for each treatment in the experiment, the matrix \code{contrasts} may have any number of columns and these are not required to be linearly independent.
Methods of assessing differential expression, such as \code{eBayes} or \code{classifyTestsF}, can then be applied to fitted model object.

The \code{coefficients} argument provides a simpler way to specify the \code{contrasts} matrix when the desired contrasts are just a subset of the original coefficients.
}
\note{
For efficiency reasons, this function does not re-factorize the design matrix for each probe.
A consequence is that, if the design matrix is non-orthogonal and the original fit included precision weights or missing values, then the unscaled standard deviations produced by this function are approximate rather than exact.
The approximation is usually acceptable.
If not, then the issue can be avoided by redefining the design matrix to fit the contrasts directly.

Even with precision weights, the results from \code{contrasts.fit} are always exact if the coefficients being compared are statistically independent.
This will always be true, for example, if the original fit was a oneway model and the group-means (no-intercept) parametrization was used for the design matrix.
}
\seealso{
An overview of linear model functions in limma is given by \link{06.LinearModels}.
}
\author{Gordon Smyth}
\examples{
#  Simulate gene expression data: 6 microarrays and 100 genes
#  with one gene differentially expressed in first 3 arrays
M <- matrix(rnorm(100*6,sd=0.3),100,6)
M[1,1:3] <- M[1,1:3] + 2
#  Design matrix corresponds to oneway layout, columns are orthogonal
design <- cbind(First3Arrays=c(1,1,1,0,0,0),Last3Arrays=c(0,0,0,1,1,1))
fit <- lmFit(M,design=design)
#  Would like to consider original two estimates plus difference between first 3 and last 3 arrays
contrast.matrix <- cbind(First3=c(1,0),Last3=c(0,1),"Last3-First3"=c(-1,1))
fit2 <- contrasts.fit(fit,contrast.matrix)
fit2 <- eBayes(fit2)
#  Large values of eb$t indicate differential expression
results <- decideTests(fit2, method="nestedF")
vennCounts(results)
}
\keyword{htest}
