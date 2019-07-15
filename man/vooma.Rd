\name{vooma}
\alias{vooma}
\alias{voomaByGroup}

\title{Convert Mean-Variance Trend to Observation-specific Precision Weights for Microarray Data}

\description{
Estimate the mean-variance relationship and use this to compute appropriate observational-level weights.
}

\usage{
vooma(y, design=NULL, correlation, block=NULL, plot=FALSE, span=NULL)
voomaByGroup(y, group, design=NULL, correlation, block=NULL,
             plot=FALSE, span=NULL, col=NULL, lwd=1, alpha=0.5,
             pch=16, cex=0.3, legend="topright")
}

\arguments{
  \item{y}{a numeric \code{matrix}, \code{EList} object, or any object containing log-expression data that can be coerced to a matrix.}
  \item{design}{design matrix with rows corresponding to samples and columns to coefficients to be estimated.  Defaults to the unit vector meaning that samples are treated as replicates.}
  \item{block}{vector or factor specifying a blocking variable on the arrays. Has length equal to the number of arrays.}
  \item{correlation}{intra-block correlation}
  \item{span}{width of the smoothing window, as a proportion of the data set.}
  \item{plot}{\code{logical} value indicating whether a plot of mean-variance trend should be displayed.}
  \item{group}{categorical vector or factor giving group membership of columns of \code{y}.}
  \item{col}{vector of colors for plotting group trends}
  \item{lwd}{line width for plotting group trends}
  \item{pch}{plotting character. Default is integer code 16 which gives a solid circle. If a vector, then should be of length \code{nrow(y)}.}
  \item{cex}{numeric vector of plot symbol expansions.  If a vector, then should be of length equal to number of groups.}
  \item{alpha}{transparancy of points, on scale from \code{0} for fully transparant to \code{1} for fully opaque.}
  \item{legend}{character string giving position to place legend.}
}

\details{
\code{vooma} is an acronym for mean-variance modelling at the observational level for arrays.

\code{vooma} estimates the mean-variance relationship in the data, and uses this to compute appropriate weights for each observation.
This done by estimating a mean-variance trend, then interpolating this trend to obtain a precision weight (inverse variance) for each observation.
The weights can then used by other functions such as \code{lmFit} to adjust for heteroscedasticity.

\code{voomaByGroup} estimates precision weights separately for each group.  In other words, it allows for different mean-variance curves in different groups.
}

\value{
An EList object with the following components:
\item{E}{numeric matrix of as input}
\item{weights}{numeric matrix of weights}
\item{design}{numeric matrix of experimental design}
\item{genes}{dataframe of gene annotation, only if \code{counts} was a \code{DGEList} object}
 }

\seealso{
\code{\link{voom}}, \code{\link{arrayWeights}}
}

\author{Charity Law and Gordon Smyth}

\references{
Law, C. (2013).
\emph{Precision weights for gene expression analysis}.
PhD Thesis. University of Melbourne, Australia.
\url{http://repository.unimelb.edu.au/10187/17598}
}
