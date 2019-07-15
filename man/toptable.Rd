\name{toptable}
\alias{toptable}
\alias{topTable}
\alias{topTableF}
\alias{topTreat}
\title{Table of Top Genes from Linear Model Fit}
\description{
Extract a table of the top-ranked genes from a linear model fit.
}
\usage{
topTable(fit, coef=NULL, number=10, genelist=fit$genes, adjust.method="BH",
         sort.by="B", resort.by=NULL, p.value=1, lfc=0, confint=FALSE)
toptable(fit, coef=1, number=10, genelist=NULL, A=NULL, eb=NULL, adjust.method="BH",
         sort.by="B", resort.by=NULL, p.value=1, lfc=0, confint=FALSE, \dots)
topTableF(fit, number=10, genelist=fit$genes, adjust.method="BH",
         sort.by="F", p.value=1, lfc=0)
topTreat(fit, coef=1, sort.by="p", resort.by=NULL, \dots)
}
\arguments{
  \item{fit}{list containing a linear model fit produced by \code{lmFit}, \code{lm.series}, \code{gls.series} or \code{mrlm}.
  For \code{topTable}, \code{fit} should be an object of class \code{MArrayLM} as produced by \code{lmFit} and \code{eBayes}.}
  \item{coef}{column number or column name specifying which coefficient or contrast of the linear model is of interest. For \code{topTable}, can also be a vector of column subscripts, in which case the gene ranking is by F-statistic for that set of contrasts.}
  \item{number}{maximum number of genes to list}
  \item{genelist}{data frame or character vector containing gene information.
  For \code{topTable} only, this defaults to \code{fit$genes}.}
  \item{A}{matrix of A-values or vector of average A-values.
  For \code{topTable} only, this defaults to \code{fit$Amean}.}
  \item{eb}{output list from \code{ebayes(fit)}.
  If \code{NULL}, this will be automatically generated.}
  \item{adjust.method}{method used to adjust the p-values for multiple testing.  Options, in increasing conservatism, include \code{"none"}, \code{"BH"}, \code{"BY"} and \code{"holm"}.
  See \code{\link{p.adjust}} for the complete list of options. A \code{NULL} value will result in the default adjustment method, which is \code{"BH"}.}
  \item{sort.by}{character string specifying statistic to rank genes by.
  Possible values for \code{topTable} and \code{toptable} are \code{"logFC"}, \code{"AveExpr"}, \code{"t"}, \code{"P"}, \code{"p"}, \code{"B"} or \code{"none"}.
  (Permitted synonyms are \code{"M"} for \code{"logFC"}, \code{"A"} or \code{"Amean"} for \code{"AveExpr"}, \code{"T"} for \code{"t"} and \code{"p"} for \code{"P"}.)
  Possibilities for \code{topTableF} are \code{"F"} or \code{"none"}.
  Possibilities for \code{topTreat} are as for \code{topTable} except for \code{"B"}.}
  \item{resort.by}{character string specifying statistic to sort the selected genes by in the output data.frame.  Possibilities are the same as for \code{sort.by}.}
  \item{p.value}{cutoff value for adjusted p-values. Only genes with lower p-values are listed.}
  \item{lfc}{minimum absolute log2-fold-change required. \code{topTable} and \code{topTableF} include only genes with (at least one) absolute log-fold-changes greater than \code{lfc}. \code{topTreat} does not remove genes but ranks genes by evidence that their log-fold-change exceeds \code{lfc}.}
  \item{confint}{logical, should confidence 95\% intervals be output for \code{logFC}?  Alternatively, can take a numeric value between zero and one specifying the confidence level required.}
  \item{\dots}{For \code{toptable}, other arguments are passed to \code{ebayes} (if \code{eb=NULL}).  For \code{topTreat}, other arguments are passed to \code{topTable}.}
}
\value{
  A dataframe with a row for the \code{number} top genes and the following columns:
  \item{genelist}{one or more columns of probe annotation, if genelist was included as input}
  \item{logFC}{estimate of the log2-fold-change corresponding to the effect or contrast (for \code{topTableF} there may be several columns of log-fold-changes)}
  \item{CI.L}{left limit of confidence interval for \code{logFC} (if \code{confint=TRUE} or \code{confint} is numeric)}
  \item{CI.R}{right limit of confidence interval for \code{logFC} (if \code{confint=TRUE} or \code{confint} is numeric)}
  \item{AveExpr}{average log2-expression for the probe over all arrays and channels, same as \code{Amean} in the \code{MarrayLM} object}
  \item{t}{moderated t-statistic (omitted for \code{topTableF})}
  \item{F}{moderated F-statistic (omitted for \code{topTable} unless more than one coef is specified)}
  \item{P.Value}{raw p-value}
  \item{adj.P.Value}{adjusted p-value or q-value}
  \item{B}{log-odds that the gene is differentially expressed (omitted for \code{topTreat})}

If \code{fit} had unique rownames, then the row.names of the above data.frame are the same in sorted order.
Otherwise, the row.names of the data.frame indicate the row number in \code{fit}.
If \code{fit} had duplicated row names, then these are preserved in the \code{ID} column of the data.frame, or in \code{ID0} if \code{genelist} already contained an \code{ID} column.
}
\details{
\code{toptable} is an earlier interface and is retained only for backward compatibility.

These functions summarize the linear model fit object produced by \code{lmFit}, \code{lm.series}, \code{gls.series} or \code{mrlm} by selecting the top-ranked genes for any given contrast.
\code{topTable} and \code{topTableF} assume that the linear model fit has already been processed by \code{\link{eBayes}}.
\code{topTreat} assumes that the fit has been processed by \code{\link{treat}}.

The p-values for the coefficient/contrast of interest are adjusted for multiple testing by a call to \code{\link[stats]{p.adjust}}.
The \code{"BH"} method, which controls the expected false discovery rate (FDR) below the specified value, is the default adjustment method because it is the most likely to be appropriate for microarray studies.
Note that the adjusted p-values from this method are bounds on the FDR rather than p-values in the usual sense.
Because they relate to FDRs rather than rejection probabilities, they are sometimes called q-values.
See \code{help("p.adjust")} for more information.

Note, if there is no good evidence for differential expression in the experiment, that it is quite possible for all the adjusted p-values to be large, even for all of them to be equal to one.
It is quite possible for all the adjusted p-values to be equal to one if the smallest p-value is no smaller than \code{1/ngenes} where \code{ngenes} is the number of genes with non-missing p-values.

The \code{sort.by} argument specifies the criterion used to select the top genes.
The choices are: \code{"logFC"} to sort by the (absolute) coefficient representing the log-fold-change; \code{"A"} to sort by average expression level (over all arrays) in descending order; \code{"T"} or \code{"t"} for absolute t-statistic; \code{"P"} or \code{"p"} for p-values; or \code{"B"} for the \code{lods} or B-statistic.

Normally the genes appear in order of selection in the output table.
If a different order is wanted, then the \code{resort.by} argument may be useful.
For example, \code{topTable(fit, sort.by="B", resort.by="logFC")} selects the top genes according to log-odds of differential expression and then orders the selected genes by log-ratio in decreasing order.
Or \code{topTable(fit, sort.by="logFC", resort.by="logFC")} would select the genes by absolute log-fold-change and then sort them from most positive to most negative.

\code{topTableF} ranks genes on the basis of moderated F-statistics for one or more coefficients.
If \code{topTable} is called and \code{coef} has two or more elements, then the specified columns will be extracted from \code{fit} and \code{topTableF} called on the result.
\code{topTable} with \code{coef=NULL} is the same as \code{topTableF}, unless the fitted model \code{fit} has only one column.

Toptable output for all probes in original (unsorted) order can be obtained by \code{topTable(fit,sort="none",n=Inf)}.
However \code{\link{write.fit}} or \code{\link{write}} may be preferable if the intention is to write the results to a file.
A related method is \code{as.data.frame(fit)} which coerces an \code{MArrayLM} object to a data.frame.

By default \code{number} probes are listed.
Alternatively, by specifying \code{p.value} and \code{number=Inf}, all genes with adjusted p-values below a specified value can be listed.

The argument \code{lfc} gives the ability to filter genes by log-fold change.
This argument is not available for \code{topTreat} because \code{treat} already handles fold-change thresholding in a more sophisticated way.
}
\note{
Although \code{topTable} enables users to set p-value and lfc cutoffs simultaneously, this is not generally recommended.
If the fold changes and p-values are not highly correlated, then the use of a fold change cutoff can increase the false discovery rate above the nominal level.
Users wanting to use fold change thresholding are usually recommended to use \code{treat} and \code{topTreat} instead.

In general, the adjusted p-values returned by \code{adjust.method="BH"} remain valid as FDR bounds only when the genes remain sorted by p-value.
Resorting the table by log-fold-change can increase the false discovery rate above the nominal level for genes at the top of resorted table.
}
\seealso{
An overview of linear model and testing functions is given in \link{06.LinearModels}.
See also \code{\link[stats]{p.adjust}} in the \code{stats} package.
}
\author{Gordon Smyth}
\examples{
#  See lmFit examples
}
\keyword{htest}
