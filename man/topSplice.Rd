\title{Top table of differentially spliced genes or exons}
\name{topSplice}
\alias{topSplice}
\description{
Top table ranking the most differentially spliced genes or exons.
}
\usage{
topSplice(fit, coef = ncol(fit), test = "simes", number = 10, FDR=1, sort.by = "p")
}
\arguments{
  \item{fit}{\code{MArrayLM} fit object produced by \code{diffSplice}.}
  \item{coef}{the coefficient (column) of fit for which differentially splicing is assessed.}
  \item{test}{character string specifying which statistical test to apply.
    Possible values are \code{"simes"}, \code{"F"} or \code{"t"}.
    \code{"F"} gives F-tests for each gene.
    \code{"t"} gives t-tests for each exon.
    \code{"simes"} gives genewise p-values derived from the t-tests after Simes adjustment for each gene.
  }
  \item{number}{integer, maximum number of rows to output.}
  \item{FDR}{numeric, only show exons or genes with false discovery rate less than this cutoff.}
  \item{sort.by}{character string specifying which column to sort results by.
    Possible values for \code{"p"}, \code{"logFC"}, \code{"NExons"} or \code{"none"}.
    \code{"logFC"} is only available if \code{test="t"} and \code{"NExons"} is only available if \code{test="simes"} or \code{test="F"}.
  }
}

\details{
Ranks genes or exons by evidence for differential splicing.
The F-statistic tests for any differences in exon usage between experimental conditions.
The exon-level t-statistics test for differences between each exon and all other exons for the same gene.

The Simes processes the exon-level p-values to give an overall call of differential splicing for each gene.
It returns the minimum Simes-adjusted p-values for each gene.

The F-tests are likely to be powerful for genes in which several exons are differentially splices.
The Simes p-values is likely to be more powerful when only a minority of the exons for a gene are differentially spliced.
The exon-level t-tests are not recommended for formal error rate control.
}

\value{A data.frame with any annotation columns found in \code{fit} plus the following columns
  \item{logFC}{log2-fold change of exon vs other exons for the same gene (if \code{level="exon"})}
  \item{t}{moderated t-statistic (if \code{level="exon"})}
  \item{F}{moderated F-statistic (if \code{level="gene"})}
  \item{P.Value}{p-value}
  \item{FDR}{false discovery rate}
}

\author{Gordon Smyth}
\seealso{
\code{\link{diffSplice}}, \code{\link{plotSplice}}

A summary of functions available in LIMMA for RNA-seq analysis is given in \link{11.RNAseq}.
}
\examples{# See diffSplice}

\keyword{rna-seq}
