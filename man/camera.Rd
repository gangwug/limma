\name{camera}
\alias{camera}
\alias{camera.default}
\alias{interGeneCorrelation}
\alias{cameraPR}
\alias{cameraPR.default}
\title{Competitive Gene Set Test Accounting for Inter-gene Correlation}
\description{
Test whether a set of genes is highly ranked relative to other genes in terms of differential expression, accounting for inter-gene correlation.
}
\usage{
\method{camera}{default}(y, index, design, contrast = ncol(design), weights = NULL,
       use.ranks = FALSE, allow.neg.cor=FALSE, inter.gene.cor=0.01, trend.var = FALSE,
       sort = TRUE, \dots)
\method{cameraPR}{default}(statistic, index, use.ranks = FALSE, inter.gene.cor=0.01, sort = TRUE, \dots)
interGeneCorrelation(y, design)
}
\arguments{
  \item{y}{a numeric matrix of log-expression values or log-ratios of expression values, or any data object containing such a matrix.
  Rows correspond to probes and columns to samples.
  Any type of object that can be processed by \code{\link{getEAWP}} is acceptable.}
  \item{statistic}{a numeric vector of genewise statistics. If \code{index} contains gene IDs, then \code{statistic} should be a named vector with the gene IDs as names.}
  \item{index}{an index vector or a list of index vectors.  Can be any vector such that \code{y[index,]} of \code{statistic[index]} selects the rows corresponding to the test set.  The list can be made using \code{\link{ids2indices}}.}
  \item{design}{design matrix.}
  \item{contrast}{contrast of the linear model coefficients for which the test is required. Can be an integer specifying a column of \code{design}, or else a numeric vector of same length as the number of columns of \code{design}.}
  \item{weights}{numeric matrix of precision weights. Can be a matrix of the same size as \code{y}, or a numeric vector of array weights with length equal to \code{ncol(y)}, or a numeric vector of gene weights with length equal to \code{nrow(y)}.}
  \item{use.ranks}{do a rank-based test (\code{TRUE}) or a parametric test (\code{FALSE})?}
  \item{allow.neg.cor}{should reduced variance inflation factors be allowed for negative correlations?}
  \item{inter.gene.cor}{numeric, optional preset value for the inter-gene correlation within tested sets.  If \code{NA} or \code{NULL}, then an inter-gene correlation will be estimated for each tested set.}
  \item{trend.var}{logical, should an empirical Bayes trend be estimated?  See \code{\link{eBayes}} for details.}
  \item{sort}{logical, should the results be sorted by p-value?}
  \item{\dots}{other arguments are not currently used}
}

\details{
\code{camera} and \code{interGeneCorrelation} implement methods proposed by Wu and Smyth (2012).
\code{camera} performs a \emph{competitive} test in the sense defined by Goeman and Buhlmann (2007).
It tests whether the genes in the set are highly ranked in terms of differential expression relative to genes not in the set.
It has similar aims to \code{geneSetTest} but accounts for inter-gene correlation.
See \code{\link{roast}} for an analogous \emph{self-contained} gene set test.

The function can be used for any microarray experiment which can be represented by a linear model.
The design matrix for the experiment is specified as for the \code{\link{lmFit}} function, and the contrast of interest is specified as for the \code{\link{contrasts.fit}} function.
This allows users to focus on differential expression for any coefficient or contrast in a linear model by giving the vector of test statistic values.

\code{camera} estimates p-values after adjusting the variance of test statistics by an estimated variance inflation factor.
The inflation factor depends on estimated genewise correlation and the number of genes in the gene set.

By default, \code{camera} uses \code{interGeneCorrelation} to estimate the mean pair-wise correlation within each set of genes.
\code{camera} can alternatively be used with a preset correlation specified by \code{inter.gene.cor} that is shared by all sets.
This usually works best with a small value, say \code{inter.gene.cor=0.01}.

If \code{interGeneCorrelation=NA}, then \code{camera} will estimate the inter-gene correlation for each set.
In this mode, \code{camera} gives rigorous error rate control for all sample sizes and all gene sets.
However, in this mode, highly co-regulated gene sets that are biological interpretable may not always be ranked at the top of the list.

With \code{interGeneCorrelation=0.01}, \code{camera} will rank biologically interpetable sets more highly.
This gives a useful compromise between strict error rate control and interpretable gene set rankings.

\code{cameraPR} is a "pre-ranked" version of \code{camera} where the genes are pre-ranked according to a pre-computed statistic.
}

\note{The default settings for \code{inter.gene.cor} and \code{allow.neg.cor} were changed to the current values in limma 3.29.6.
Previously, the default was to estimate an inter-gene correlation for each set.
To reproduce the earlier default, use \code{allow.neg.cor=TRUE} and \code{inter.gene.cor=NA}.}

\value{
\code{camera} and \code{cameraPR} return a data.frame with a row for each set and the following columns:
\item{NGenes}{number of genes in set.}
\item{Correlation}{inter-gene correlation (only included if the \code{inter.gene.cor} was not preset).}
\item{Direction}{direction of change (\code{"Up"} or \code{"Down"}).}
\item{PValue}{two-tailed p-value.}
\item{FDR}{Benjamini and Hochberg FDR adjusted p-value.}

\code{interGeneCorrelation} returns a list with components:
\item{vif}{variance inflation factor.}
\item{correlation}{inter-gene correlation.}
}

\author{Di Wu and Gordon Smyth}

\references{
Wu, D, and Smyth, GK (2012). Camera: a competitive gene set test accounting for inter-gene correlation.
\emph{Nucleic Acids Research} 40, e133.
\url{http://nar.oxfordjournals.org/content/40/17/e133}

Goeman, JJ, and Buhlmann, P (2007).
Analyzing gene expression data in terms of gene sets: methodological issues.
\emph{Bioinformatics} 23, 980-987. 
}

\seealso{
\code{\link{getEAWP}}

\code{\link{rankSumTestWithCorrelation}},
\code{\link{geneSetTest}},
\code{\link{roast}},
\code{\link{fry}},
\code{\link{romer}},
\code{\link{ids2indices}}.

There is a topic page on \link{10.GeneSetTests}.
}

\examples{
y <- matrix(rnorm(1000*6),1000,6)
design <- cbind(Intercept=1,Group=c(0,0,0,1,1,1))

# First set of 20 genes are genuinely differentially expressed
index1 <- 1:20
y[index1,4:6] <- y[index1,4:6]+1

# Second set of 20 genes are not DE
index2 <- 21:40
 
camera(y, index1, design)
camera(y, index2, design)

camera(y, list(set1=index1,set2=index2), design, inter.gene.cor=NA)
camera(y, list(set1=index1,set2=index2), design, inter.gene.cor=0.01)

# Pre-ranked version
fit <- eBayes(lmFit(y, design))
cameraPR(fit$t[,2], list(set1=index1,set2=index2))
}

\keyword{gene set test}
\concept{gene set test}
\concept{gene set enrichment analysis}
