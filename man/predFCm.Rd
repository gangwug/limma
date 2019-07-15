\name{predFCm}
\alias{predFCm}
\alias{predFCm}
\title{Predictive log fold change for microarrays}
\description{
Calculate the predictive log fold change for a particular coefficient from a fit object.
}
\usage{
predFCm(fit, coef=2, var.indep.of.fc=TRUE, all.de=TRUE, prop.true.null.method="lfdr")
}

\arguments{
 \item{fit}{an \code{MArrayLM} fitted model object produced by \code{lmFit} and \code{eBayes}}
 \item{coef}{integer vector indicating which columns in the fit object are to be shrunk}
 \item{var.indep.of.fc}{assume the genewise variances are independent of genewise fold changes?}
 \item{all.de}{assume all genes are have a non-zero true fold change (\code{TRUE})?
 If \code{FALSE}, then the proportion of truly non-differentially (non-DE) genes expressed will be estimated.}
 \item{prop.true.null.method}{method used to estimate proportion of truly non-DE genes. See \code{\link{propTrueNull}} for possible values.}
}

\details{
The predictive log fold changes are calculated as the posterior mean log fold changes in the empirical Bayes hierarchical model.
We call them predictive log fold changes because they are the best prediction of what the log fold change will be for each gene in a comparable future experiment.

The log fold changes are shrunk towards zero depending on how variable they are.
The \code{var.indep.of.fc} argument specifies whether the prior belief is that the log fold changes are independent of the variability of the genes or whether the log fold changes increase with increasing variability of the genes.

If \code{all.de=TRUE}, then all genes are assumed to have a non-zero log fold change, even if quite small.
If \code{all.de=FALSE}, then some genes are assumed to have log fold changes exactly zero.
The proportion of non-DE genes is estimated and taken into account in the calculation.
}

\value{numeric vector of predictive (shrunk) log fold changes}

\seealso{
\code{\link{lmFit}}, \code{\link{eBayes}}, \code{\link{contrasts.fit}}
}

\author{Belinda Phipson and Gordon Smyth}

\references{
Phipson, B. (2013).
\emph{Empirical Bayes modelling of expression profiles and their associations}.
PhD Thesis. University of Melbourne, Australia.
\url{http://repository.unimelb.edu.au/10187/17614}
}

\examples{
#  Simulate gene expression data,
#  6 microarrays with 1000 genes on each array 
set.seed(2004)
y <- matrix(rnorm(6000),ncol=4)

# two experimental groups and one control group with two replicates each
group <- factor(c("A","A","B","B"))
design <- model.matrix(~group)

# fit a linear model
fit <- lmFit(y,design)
fit <- eBayes(fit)

# output predictive log fold changes for first 5 genes
pfc <- predFCm(fit,coef=2)
}
