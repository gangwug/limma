\title{Differential splicing plot with junctions}
\name{plotExonJunc}
\alias{plotExonJunc}
\description{
Plot differential usage results by exons and junctions for the specified gene and highlight the significantly spliced exons and junctions.
}
\usage{
plotExonJunc(fit, coef=ncol(fit), geneid, genecolname=NULL, FDR=0.05, annotation=NULL)
}
\arguments{
  \item{fit}{\code{MArrayLM} fit object produced by \code{diffSplice}. Must have the Entrez gene ids for all the exons and junctions stored in \code{fit$genes$GeneID}, length information for all the exons and junctions stored in \code{fit$genes$Length} and the strand information stored in \code{fit$genes$Strand}. To distinguish between exons and junctions, \code{fit$genes$Length} has to be set to 1 for all the junctions.}
  \item{coef}{the coefficient (column) of fit for which differentially splicing is assessed.}
  \item{geneid}{character string, ID of the gene to plot.}
  \item{genecolname}{column name of \code{fit$genes} containing \code{geneid}.}
  \item{FDR}{numeric, highlight exons and junctions with false discovery rate less than this cutoff. Red indicates up-regulation whereas blue indicates down-regulation. The FDR of the individual exon/junction is calculated based on the exon-level t-statistics test for differences between each exon/junction and all other exons/junctions for the same gene.}
  \item{annotation}{data frame containing the full exon annotation of the corresponding species. Must have the Entrez gene ids for all the exons stored in the \code{GeneID} column, start and end positions for all the exons stored in the \code{Start} and \code{End} columns, respectively.}
}

\details{
Plot differential usage results by exons and junctions for the specified gene.
The significantly spliced individual exons are highlighted as red blocks if up-regulated and blue blocks if down-regulated. 
All other exons are displayed as black blocks.
The significantly spliced individual junctions are highlighted as red lines if up-regulated and blue lines if down-regulated.
All other junctions are displayed as black lines.

Since the \code{diffSplice} analysis is usually performed after filtering, the full annotation (e.g. the inbuilt annotation in \code{featureCounts}) is highly recommended for producing the plot. When \code{annotation} is provided, the filtered exons are displayed as grey blocks.
}

\value{A plot is created on the current graphics device.}
\author{Yunshun Chen and Gordon Smyth}
\seealso{
\code{\link{diffSplice}}, \code{\link{topSplice}}
}
\examples{
\dontrun{
# diffSplice analysis
v <- voom(dge, design)
fit <- lmFit(v, design)
ex <- diffSplice(fit, geneid="GeneID")

# Get full annotation from Rsubread
library(Rsubread)
annotation.full <- getInBuiltAnnotation("mm10")

# Make a plot
plotExonJunc(ex, geneid="Foxp1", genecolname="Symbol", annotation=annotation.full)
}
}

\keyword{hplot}
\keyword{rna-seq}
