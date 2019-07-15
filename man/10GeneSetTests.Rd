\name{10.GeneSetTests}
\alias{10.GeneSetTests}
\title{Topic: Gene Set Tests}

\description{
This page gives an overview of the LIMMA functions for gene set testing and pathway analysis.

\describe{
\item{ \code{\link{roast}} }{
	Self-contained gene set testing for one set.}

\item{ \code{\link{mroast}} }{
	Self-contained gene set testing for many sets.}

\item{ \code{\link{fry}} }{
	Fast approximation to \code{mroast}, especially useful when heteroscedasticity of genes can be ignored.}

\item{ \code{\link{camera}} }{
	Competitive gene set testing.}

\item{ \code{\link{romer}} and \code{\link{topRomer}} }{
	Gene set enrichment analysis.}

\item{ \code{\link{ids2indices}} }{
	Convert gene sets consisting of vectors of gene identifiers into a list of indices suitable for use in the above functions.}

\item{ \code{\link{alias2Symbol}} and \code{\link{alias2SymbolTable}} }{
	Convert gene symbols or aliases to current official symbols.}

\item{ \code{\link{geneSetTest}} or \code{\link{wilcoxGST}} }{
	Simple gene set testing based on gene or probe permutation.}

\item{ \code{\link{barcodeplot}} }{
	Enrichment plot of a gene set.}

\item{ \code{\link{goana}} and \code{\link{topGO}}}{
	Gene ontology over-representation analysis of gene lists using Entrez Gene IDs.
	\code{goana} can work directly on a fitted model object or on one or more lists of genes.}

\item{ \code{\link{kegga}} and \code{\link{topKEGG}}}{
	KEGG pathway over-representation analysis of gene lists using Entrez Gene IDs.
	\code{kegga} can work directly on a fitted model object or on one or more lists of genes.}
}
}

\author{Gordon Smyth}

\seealso{
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
