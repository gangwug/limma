\name{11.RNAseq}
\alias{11.RNAseq}
\title{Topic: Analysis of RNA-seq Data}

\description{
This page gives an overview of LIMMA functions to analyze RNA-seq data.

\describe{
\item{ \code{\link{voom}} }{
	Transform RNA-seq or ChIP-seq counts to log counts per million (log-cpm) with associated precision weights.
	After this tranformation, RNA-seq or ChIP-seq data can be analyzed using the same functions as would be used for microarray data.}

\item{ \code{\link{voomWithQualityWeights}} }{
	Combines the functionality of \code{voom} and \code{arrayWeights}.}

\item{ \code{\link{diffSplice}} }{
	Test for differential exon usage between experimental conditions.}

\item{ \code{\link{topSplice}} }{
	Show a data.frame of top results from \code{diffSplice}.}

\item{ \code{\link{plotSplice}} }{
	Plot results from \code{diffSplice}.}

\item{ \code{\link{plotExons}} }{
	Plot logFC for individual exons for a given gene.}
}
}

\references{
Law, CW, Chen, Y, Shi, W, Smyth, GK (2014).
Voom: precision weights unlock linear model analysis tools for RNA-seq read counts.
\emph{Genome Biology} 15, R29.
\url{http://genomebiology.com/2014/15/2/R29}

Ritchie, ME, Phipson, B, Wu, D, Hu, Y, Law, CW, Shi, W, and Smyth, GK (2015).
limma powers differential expression analyses for RNA-sequencing and microarray studies.
\emph{Nucleic Acids Research} 43, e47.
\url{http://nar.oxfordjournals.org/content/43/7/e47}
}

\seealso{
See also the edgeR package for normalization and data summaries of RNA-seq data, as well as for alternative differential expression methods based on the negative binomial distribution.
\code{voom} accepts DGEList objects and normalization factors from edgeR.

\link{01.Introduction},
\link{02.Classes},
\link{03.ReadingData},
\link{04.Background},
\link{05.Normalization},
\link{06.LinearModels},
\link{07.SingleChannel},
\link{08.Tests},
\link{09.Diagnostics},
\link{10.GeneSetTests},
\link{11.RNAseq}
}

\keyword{documentation}
