\name{diffSplice}
\alias{diffSplice}
\title{Test for Differential Splicing}
\description{Given a linear model fit at the exon level, test for differences in exon retention between experimental conditions.}
\usage{
diffSplice(fit, geneid, exonid=NULL, robust=FALSE, verbose=TRUE)
}
\arguments{
  \item{fit}{an \code{MArrayLM} fitted model object produced by \code{lmFit} or \code{contrasts.fit}. Rows should correspond to exons.}
  \item{geneid}{gene identifiers. Either a vector of length \code{nrow(fit)} or the name of the column of \code{fit$genes} containing the gene identifiers. Rows with the same ID are assumed to belong to the same gene.}
  \item{exonid}{exon identifiers. Either a vector of length \code{nrow(fit)} or the name of the column of \code{fit$genes} containing the exon identifiers.}
  \item{robust}{logical, should the estimation of the empirical Bayes prior parameters be robustified against outlier sample variances?}
  \item{verbose}{logical, if \code{TRUE} some diagnostic information about the number of genes and exons is output.}
}
\value{
An object of class \code{MArrayLM} containing both exon level and gene level tests.
Results are sorted by geneid and by exonid within gene.
  \item{coefficients}{numeric matrix of coefficients of same dimensions as \code{fit}. Each coefficient is the difference between the log-fold-change for that exon versus the average log-fold-change for all other exons for the same gene.}
  \item{t}{numeric matrix of moderated t-statistics, of same dimensions as \code{fit}.}
  \item{p.value}{numeric vector of p-values corresponding to the t-statistics}
  \item{genes}{data.frame of exon annotation}
  \item{genecolname}{character string giving the name of the column of \code{genes} containing gene IDs}
  \item{gene.F}{numeric matrix of moderated F-statistics, one row for each gene.}
  \item{gene.F.p.value}{numeric matrix of p-values corresponding to \code{gene.F}}
  \item{gene.simes.p.value}{numeric matrix of Simes adjusted p-values, one row for each gene.}
  \item{gene.genes}{data.frame of gene annotation.}
}
\details{
This function tests for differential exon usage for each gene and for each column of \code{fit}.

Testing for differential exon usage is equivalent to testing whether the log-fold-changes in the \code{fit} differ between exons for the same gene.
Two different tests are provided.
The first is an F-test for differences between the log-fold-changes.
The other is a series of t-tests in which each exon is compared to the average of all other exons for the same gene.
The exon-level t-tests are converted into a genewise test by adjusting the p-values for the same gene by Simes method.
The minimum adjusted p-value is then used for each gene.

This function can be used on data from an exon microarray or can be used in conjunction with voom for exon-level RNA-seq counts.
}
\seealso{
\code{\link{topSplice}}, \code{\link{plotSplice}}

A summary of functions available in LIMMA for RNA-seq analysis is given in \link{11.RNAseq}.
}
\author{Gordon Smyth and Charity Law}

\examples{
\dontrun{
v <- voom(dge,design)
fit <- lmFit(v,design)
ex <- diffSplice(fit,geneid="EntrezID")
topSplice(ex)
plotSplice(ex)
}
}

\keyword{rna-seq}
