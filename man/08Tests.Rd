\name{08.Tests}
\alias{08.Tests}
\title{Topic: Hypothesis Testing for Linear Models}

\description{
LIMMA provides a number of functions for multiple testing across both contrasts and genes.
The starting point is an \code{MArrayLM} object, called \code{fit} say, resulting from fitting a linear model and running \code{eBayes} and, optionally, \code{contrasts.fit}.
See \link{06.LinearModels} or \link{07.SingleChannel} for details.
}

\section{Multiple testing across genes and contrasts}{
The key function is \code{\link{decideTests}}.
This function writes an object of class \code{\link[limma:TestResults]{TestResults}}, which is basically a matrix of \code{-1}, \code{0} or \code{1} elements, of the same dimension as \code{fit$coefficients}, indicating whether each coefficient is significantly different from zero.
A number of different multiple testing strategies are provided.
\code{decideTests} calls \code{\link{classifyTestsF}} to implement the nested F-test strategt. 

\code{\link{selectModel}} chooses between linear models for each probe using AIC or BIC criteria.
This is an alternative to hypothesis testing and can choose between non-nested models.

A number of other functions are provided to display the results of \code{decideTests}.
The functions \code{\link{heatDiagram}} (or the older version \code{\link{heatdiagram}} displays the results in a heat-map style display.
This allows visual comparison of the results across many different conditions in the linear model.

The functions \code{\link{vennCounts}} and \code{\link{vennDiagram}} provide Venn diagram style summaries of the results.

Summary and \code{show} method exists for objects of class \code{TestResults}.

The results from \code{decideTests} can also be included when the results of a linear model fit are written to a file using \code{\link{write.fit}}.
}

\section{Gene Set Tests}{
Competitive gene set testing for an individual gene set is provided by \code{\link{wilcoxGST}} or \code{\link{geneSetTest}}, which permute genes.
The gene set can be displayed using \code{\link{barcodeplot}}.

Self-contained gene set testing for an individual set is provided by \code{\link{roast}}, which uses rotation technology, analogous to permuting arrays.

Gene set enrichment analysis for a large database of gene sets is provided by \code{\link{romer}}.
\code{\link{topRomer}} is used to rank results from \code{romer}.

The functions \code{\link{alias2Symbol}}, \code{\link{alias2SymbolTable}} and \code{\link{alias2SymbolUsingNCBI}} are provided to help match gene sets with microarray probes by way of official gene symbols.
}

\section{Global Tests}{
The function \code{\link{genas}} can test for associations between two contrasts in a linear model.

Given a set of p-values, the function \code{\link{propTrueNull}} can be used to estimate the proportion of true null hypotheses.

When evaluating test procedures with simulated or known results, the utility function \code{\link{auROC}} can be used to compute the area under the Receiver Operating Curve for the test results for a given probe.
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
