\title{Plot exons of differentially expressed gene}
\name{plotExons}
\alias{plotExons}
\description{
Plot exons of differentially expressed gene and mark the differentially expressed exons.
}
\usage{
plotExons(fit, coef = ncol(fit), geneid = NULL, genecolname = "GeneID",
          exoncolname = NULL, rank = 1L, FDR = 0.05)
}
\arguments{
  \item{fit}{\code{MArrayLM} fit object produced by \code{eBayes}.}
  \item{coef}{the coefficient (column) of fit for which differential expression is assessed.}
  \item{geneid}{character string, ID of the gene to plot.}
  \item{genecolname}{character string for the column name of \code{fit$genes} containing gene IDs. Defaults to "GeneID" for Entrez Gene ID.}
  \item{exoncolname}{character string for the column name of \code{fit$genes} containing exon IDs.}
  \item{rank}{integer, if \code{geneid=NULL} then this ranked gene will be plotted.}
  \item{FDR}{numeric, mark differentially expressed exons with false discovery rate less than this cutoff.}
}

\details{
Plots log2-fold-change by exon for the specified gene and highlight the differentially expressed exons.
Show annotations such as GeneID, Symbol and Strand if available as title for the gene to plot.
The significantly differentially expressed individual exons are highlighted as red dots for up-regulation and as blue dots for down-regulation.
The size of the dots are weighted by its significance.
}

\value{A plot is created on the current graphics device.}
\author{Yifang Hu and Gordon Smyth}
\seealso{
\code{\link{lmFit}}, \code{\link{eBayes}}, \code{\link{plotSplice}}

A summary of functions available in LIMMA for RNA-seq analysis is given in \link{11.RNAseq}.
}
\examples{

\dontrun{
fit <- lmFit(y,design)
fit <- eBayes(fit)
plotExons(fit)
plotExons(fit, exoncolname = "Start", rank = 1)
plotExons(fit, geneid = "ps", genecolname = "Symbol", exoncolname = "Start")
}

}

\keyword{hplot}
\keyword{rna-seq}
