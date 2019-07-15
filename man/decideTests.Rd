\name{decideTests}
\alias{decideTests}
\alias{decideTests.default}
\alias{decideTests.MArrayLM}
\title{Multiple Testing Across Genes and Contrasts}
\description{
Identify which genes are significantly differentially expressed for each contrast from a fit object containing p-values and test statistics.
A number of different multiple testing strategies are offered that adjust for multiple testing down the genes as well as across contrasts for each gene.
}
\usage{
\method{decideTests}{MArrayLM}(object, method = "separate", adjust.method = "BH", p.value = 0.05,
            lfc = 0, \dots)
\method{decideTests}{default}(object, method = "separate", adjust.method = "BH", p.value = 0.05,
            lfc = 0, coefficients = NULL, cor.matrix = NULL, tstat = NULL, df = Inf,
            genewise.p.value = NULL, \dots)
}
\arguments{
  \item{object}{a numeric matrix of p-values or an \code{MArrayLM} object from which p-values and t-statistics can be extracted.}
  \item{method}{character string specifying how genes and contrasts are to be combined in the multiple testing scheme.  Choices are \code{"separate"}, \code{"global"}, \code{"hierarchical"} or \code{"nestedF"}.}
  \item{adjust.method}{character string specifying p-value adjustment method.  Possible values are \code{"none"}, \code{"BH"}, \code{"fdr"} (equivalent to \code{"BH"}), \code{"BY"} and \code{"holm"}. See \code{\link[stats]{p.adjust}} for details.}
  \item{p.value}{numeric value between 0 and 1 giving the required family-wise error rate or false discovery rate.}
  \item{lfc}{numeric, minimum absolute log2-fold-change required.}
  \item{coefficients}{numeric matrix of coefficients or log2-fold-changes. Of same dimensions as \code{object}.}
  \item{cor.matrix}{correlation matrix of coefficients. Square matrix of dimension \code{ncol(object)}.}
  \item{tstat}{numeric matrix of t-statistics. Of same dimensions as \code{object}.}
  \item{df}{numeric vector of length \code{nrow(object)} giving degrees of freedom for the t-statistics.}
  \item{genewise.p.value}{numeric vector of length \code{nrow(object)} containing summary gene-level p-values for use with \code{method="hierarchical"}.}
  \item{\dots}{other arguments are not used.}
}
\value{
An object of class \code{\link[=TestResults-class]{TestResults}}.
This is essentially a numeric matrix with elements \code{-1}, \code{0} or \code{1} depending on whether each t-statistic is classified as significantly negative, not significant or significantly positive.

If \code{lfc>0} then contrasts are judged significant only when the log2-fold change is at least this large in absolute value.
For example, one might choose \code{lfc=log2(1.5)} to restrict to 50\% changes or \code{lfc=1} for 2-fold changes.
In this case, contrasts must satisfy both the p-value and the fold-change cutoff to be judged significant.
}
\details{
This function can be applied to a matrix of p-values but is more often applied to an \code{MArrayLM} fit object produced by \code{eBayes} or \code{treat}.
In either case, rows of \code{object} correspond to genes and columns to coefficients or contrasts.

This function applies a multiple testing procedure and a significance level cutoff to the statistics contained in \code{object}.
It implements a number of multiple testing procedures for determining whether each statistic should be considered significantly different from zero.

\code{method="separate"} will apply multiple testing adjustments to each column of p-values separately.
Setting \code{method="separate"} is equivalent to using \code{topTable} separately for each coefficient in the linear model fit and will identify the same probes as significantly differentially expressed if \code{adjust.method} is the same.
\code{method="global"} will treat the entire matrix of t-statistics as a single vector of unrelated tests.
\code{method="hierarchical"} adjusts down genes and then across contrasts.
\code{method="nestedF"} adjusts down genes and then uses \code{classifyTestsF} to classify contrasts as significant or not for the selected genes.

The default \code{method="separate"} and \code{adjust.method="BH"} settings are appropriate for most analyses.
\code{method="global"} is useful when it is important that the same t-statistic cutoff should correspond to statistical significance for all the contrasts.
The \code{"nestedF"} method was proposed by Michaud et al (2008) and achieves better consistency between contrasts than the other methods.
It provides formal error rate control at the gene level but not for individual contrasts.
Please see the limma User's Guide for a discussion of the statistical properties of the various methods.
}
\note{
Although this function enables users to set p-value and lfc cutoffs simultaneously, this combination criterion is not generally recommended.
Unless the fold changes and p-values are very highly correlated, the addition of a fold change cutoff can increase the family-wise error rate or false discovery rate above the nominal level.
Users wanting to use fold change thresholding are recommended to use \code{treat} instead of \code{eBayes} and to leave \code{lfc} at the default value when using \code{decideTests}.
}
\references{
Michaud, J, Simpson, KM, Escher, R, Buchet-Poyau, K, Beissbarth, T, Carmichael, C, Ritchie, ME, Schutz, F, Cannon, P, Liu, M, Shen, X, Ito, Y, Raskind, WH, Horwitz, MS, Osato, M, Turner, DR, Speed, TP, Kavallaris, M, Smyth, GK, and Scott, HS (2008).
Integrative analysis of RUNX1 downstream pathways and target genes.
\emph{BMC Genomics} 9, 363. 
}
\seealso{
An overview of multiple testing functions is given in \link{08.Tests}.
}
\author{Gordon Smyth}
\keyword{htest}
