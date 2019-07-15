\name{fitmixture}
\alias{fitmixture}
\title{Fit Mixture Model by Non-Linear Least Squares}
\description{Fit Mixture Model by Non-Linear Least Squares}

\usage{
fitmixture(log2e, mixprop, niter = 4, trace = FALSE)
}

\arguments{
  \item{log2e}{a numeric matrix containing log2 expression values. Rows correspond to probes for genes and columns to RNA samples.}
  \item{mixprop}{a vector of length \code{ncol(log2e)} giving the mixing proportion (between 0 and 1) for each sample.}
  \item{niter}{integer number of iterations.}
  \item{trace}{logical. If \code{TRUE}, summary working estimates are output from each iteration.}
}

\details{
A mixture experiment is one in which two reference RNA sources are mixed in different proportions to create experimental samples.
Mixture experiments have been used to evaluate genomic technologies and analysis methods (Holloway et al, 2006).
This function uses all the data for each gene to estimate the expression level of the gene in each of two pure samples.

The function fits a nonlinear mixture model to the log2 expression values for each gene.
The expected values of \code{log2e} for each gene are assumed to be of the form
\code{log2( mixprop*Y1 + (1-mixprop)*Y2 )}
where \code{Y1} and \code{Y2} are the expression levels of the gene in the two reference samples being mixed.
The \code{mixprop} values are the same for each gene but \code{Y1} and \code{Y2} are specific to the gene.
The function returns the estimated values \code{A=0.5*log2(Y1*Y2)} and \code{M=log2(Y2/Y1)} for each gene.

The nonlinear estimation algorithm implemented in \code{fitmixture} uses a nested Gauss-Newton iteration (Smyth, 1996).
It is fully vectorized so that the estimation is done for all genes simultaneously.
}

\value{List with three components:
\item{A}{numeric vector giving the estimated average log2 expression of the two reference samples for each gene}
\item{M}{numeric vector giving estimated log-ratio of expression between the two reference samples for each gene}
\item{stdev}{standard deviation of the residual term in the mixture model for each gene}
}

\author{Gordon K Smyth}

\references{
Holloway, A. J., Oshlack, A., Diyagama, D. S., Bowtell, D. D. L., and Smyth, G. K. (2006).
Statistical analysis of an RNA titration series evaluates microarray precision and sensitivity on a whole-array basis.
\emph{BMC Bioinformatics} 7, Article 511.
\url{http://www.biomedcentral.com/1471-2105/7/511}

Smyth, G. K. (1996).
Partitioned algorithms for maximum likelihood and other nonlinear estimation.
\emph{Statistics and Computing}, 6, 201-216. 
\url{http://www.statsci.org/smyth/pubs/partitio.pdf}
}

\examples{
ngenes <- 100
TrueY1 <- rexp(ngenes)
TrueY2 <- rexp(ngenes)
mixprop <- matrix(c(0,0.25,0.75,1),1,4)
TrueExpr <- TrueY1 %*% mixprop + TrueY2 %*% (1-mixprop)

log2e <- log2(TrueExpr) + matrix(rnorm(ngenes*4),ngenes,4)*0.1
out <- fitmixture(log2e,mixprop)

# Plot true vs estimated log-ratios
plot(log2(TrueY1/TrueY2), out$M)
}

