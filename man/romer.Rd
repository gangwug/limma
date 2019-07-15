\name{romer}
\alias{romer}
\alias{romer.default}

\title{Rotation Gene Set Enrichment Analysis}
\description{
Gene set enrichment analysis for linear models using rotation tests (ROtation testing using MEan Ranks).
}
\usage{
\method{romer}{default}(y, index, design = NULL, contrast = ncol(design),
      array.weights = NULL, block = NULL, correlation,
      set.statistic = "mean", nrot = 9999, shrink.resid = TRUE, \dots)
}
\arguments{
  \item{y}{numeric matrix giving log-expression values.}
  \item{index}{list of indices specifying the rows of \code{y} in the gene sets. The list can be made using \code{\link{ids2indices}}.}
  \item{design}{design matrix.}
  \item{contrast}{contrast for which the test is required. Can be an integer specifying a column of \code{design}, or else a contrast vector of length equal to the number of columns of \code{design}.}
  \item{array.weights}{optional numeric vector of array weights.} 
  \item{block}{optional vector of blocks.}
  \item{correlation}{correlation between blocks.}
  \item{set.statistic}{statistic used to summarize the gene ranks for each set. Possible values are \code{"mean"}, \code{"floormean"} or \code{"mean50"}.}
  \item{nrot}{number of rotations used to estimate the p-values.}
  \item{shrink.resid}{logical, should the residuals be shrunk to remove systematics effects before rotation.}
  \item{\dots}{other arguments not currently used.}
}
\value{
Numeric matrix giving p-values and the number of matched genes in each gene set.
Rows correspond to gene sets.
There are four columns giving the number of genes in the set and p-values for the alternative hypotheses mixed, up or down.
}
\details{
This function implements the ROMER procedure described by Majewski et al (2010) and Ritchie et al (2015).
\code{romer} tests a hypothesis similar to that of Gene Set Enrichment Analysis (GSEA) (Subramanian et al, 2005) but is designed for use with linear models.
Like GSEA, it is designed for use with a database of gene sets.
Like GSEA, it is a competitive test in that the different gene sets are pitted against one another.
Instead of permutation, it uses rotation, a parametric resampling method suitable for linear models (Langsrud, 2005; Wu et al, 2010).
\code{romer} can be used with any linear model with some level of replication.

In the output, p-values are given for each set for three possible alternative hypotheses.
The alternative "up" means the genes in the set tend to be up-regulated, with positive t-statistics.
The alternative "down" means the genes in the set tend to be down-regulated, with negative t-statistics.
The alternative "mixed" test whether the genes in the set tend to be differentially expressed, without regard for direction.
In this case, the test will be significant if the set contains mostly large test statistics, even if some are positive and some are negative.
The first two alternatives are appropriate if you have a prior expection that all the genes in the set will react in the same direction.
The "mixed" alternative is appropriate if you know only that the genes are involved in the relevant pathways, without knowing the direction of effect for each gene.

Note that \code{romer} estimates p-values by simulation, specifically by random rotations of the orthogonalized residuals (called effects in R).
This means that the p-values will vary slightly from run to run.
To get more precise p-values, increase the number of rotations \code{nrot}.
By default, the orthogonalized residual corresponding to the contrast being tested is shrunk have the same expected squared size as a null residual.

The argument \code{set.statistic} controls the way that t-statistics are summarized to form a summary test statistic for each set.
In all cases, genes are ranked by moderated t-statistic.
If \code{set.statistic="mean"}, the mean-rank of the genes in each set is the summary statistic.
If \code{set.statistic="floormean"} then negative t-statistics are put to zero before ranking for the up test, and vice versa for the down test.
This improves the power for detecting genes with a subset of responding genes.
If \code{set.statistics="mean50"}, the mean of the top 50\% ranks in each set is the summary statistic.
This statistic performs well in practice but is slightly slower to compute.
See Wu et al (2010) for discussion of these set statistics.
}

\seealso{
\code{\link{topRomer}},
\code{\link{ids2indices}},
\code{\link{roast}},
\code{\link{camera}},
\code{\link{wilcoxGST}}

There is a topic page on \link{10.GeneSetTests}.
}
\author{Yifang Hu and Gordon Smyth}


\references{
Langsrud, O (2005).
Rotation tests.
\emph{Statistics and Computing} 15, 53-60

Majewski, IJ, Ritchie, ME, Phipson, B, Corbin, J, Pakusch, M, Ebert, A, Busslinger, M, Koseki, H, Hu, Y, Smyth, GK, Alexander, WS, Hilton, DJ, and Blewitt, ME (2010).
Opposing roles of polycomb repressive complexes in hematopoietic stem and progenitor cells.
\emph{Blood} 116, 731-739.
\url{http://www.ncbi.nlm.nih.gov/pubmed/20445021}

Ritchie, ME, Phipson, B, Wu, D, Hu, Y, Law, CW, Shi, W, and Smyth, GK (2015).
limma powers differential expression analyses for RNA-sequencing and microarray studies.
\emph{Nucleic Acids Research} 43, e47.
\url{http://nar.oxfordjournals.org/content/43/7/e47}

Subramanian, A, Tamayo, P, Mootha, VK, Mukherjee, S, Ebert, BL, Gillette, MA, Paulovich, A, Pomeroy, SL, Golub, TR, Lander, ES and Mesirov JP (2005).
Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles.
\emph{Proc Natl Acad Sci U S A} 102, 15545-15550

Wu, D, Lim, E, Francois Vaillant, F, Asselin-Labat, M-L, Visvader, JE, and Smyth, GK (2010).
ROAST: rotation gene set tests for complex microarray experiments.
\emph{Bioinformatics} 26, 2176-2182.
\url{http://bioinformatics.oxfordjournals.org/content/26/17/2176}
}

\examples{
y <- matrix(rnorm(100*4),100,4)
design <- cbind(Intercept=1,Group=c(0,0,1,1))
index <- 1:5
y[index,3:4] <- y[index,3:4]+3

index1 <- 1:5
index2 <- 6:10
r <- romer(y=y,index=list(set1=index1,set2=index2),design=design,contrast=2,nrot=99)
r
topRomer(r,alt="up")
topRomer(r,alt="down")
}

\keyword{gene set test}
\concept{gene set test}
\concept{gene set enrichment analysis}
