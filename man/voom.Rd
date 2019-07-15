\name{voom}
\alias{voom}
\title{Transform RNA-Seq Data Ready for Linear Modelling}
\description{
Transform count data to log2-counts per million (logCPM), estimate the mean-variance relationship and use this to compute appropriate observation-level weights.
The data are then ready for linear modelling.
}

\usage{
voom(counts, design = NULL, lib.size = NULL, normalize.method = "none",
     span = 0.5, plot = FALSE, save.plot = FALSE, \dots)
}
\arguments{
 \item{counts}{a numeric \code{matrix} containing raw counts, or an \code{ExpressionSet} containing raw counts, or a \code{DGEList} object. Counts must be non-negative and NAs are not permitted.}
 \item{design}{design matrix with rows corresponding to samples and columns to coefficients to be estimated.  Defaults to the unit vector meaning that samples are treated as replicates.}
 \item{lib.size}{numeric vector containing total library sizes for each sample.
 Defaults to the normalized (effective) library sizes in \code{counts} if \code{counts} is a \code{DGEList} or
 to the columnwise count totals if \code{counts} is a matrix.}
 \item{normalize.method}{the microarray-style normalization method to be applied to the logCPM values (if any).
 Choices are as for the \code{method} argument of \code{normalizeBetweenArrays} when the data is single-channel.
 Any normalization factors found in \code{counts} will still be used even if \code{normalize.method="none"}.}
 \item{span}{width of the lowess smoothing window as a proportion.}
 \item{plot}{logical, should a plot of the mean-variance trend be displayed?}
 \item{save.plot}{logical, should the coordinates and line of the plot be saved in the output?}
 \item{\dots}{other arguments are passed to \code{lmFit}.}
  }

\details{
This function is intended to process RNA-Seq or ChIP-Seq data prior to linear modelling in limma.

\code{voom} is an acronym for mean-variance modelling at the observational level.
The idea is to estimate the mean-variance relationship in the data, then use this to compute an appropriate precision weight for each observation.
Count data always show marked mean-variance relationships.
Raw counts show increasing variance with increasing count size, while log-counts typically show a decreasing mean-variance trend.
This function estimates the mean-variance trend for log-counts, then assigns a weight to each observation based on its predicted variance.
The weights are then used in the linear modelling process to adjust for heteroscedasticity. 

\code{voom} performs the following specific calculations.
First, the counts are converted to logCPM values, adding 0.5 to all the counts to avoid taking the logarithm of zero.
The matrix of logCPM values is then optionally normalized.
The \code{lmFit} function is used to fit row-wise linear models.
The \code{lowess} function is then used to fit a trend to the square-root-standard-deviations as a function of an average log-count measure.
The trend line is then used to predict the variance of each logCPM value as a function of its fitted value on the count scale, and the inverse variances become the estimated precision weights.

For good results, the \code{counts} matrix should be filtered to remove remove rows with very low counts before running voom().
The \code{filterByExpr} function in the edgeR package can be used for that purpose.

If \code{counts} is a \code{DGEList} object from the edgeR package, then voom will use the normalization factors found in the object when computing the logCPM values.
In other words, the logCPM values are computed from the effective library sizes rather than the raw library sizes.
If the \code{DGEList} object has been scale-normalized in edgeR, then it is usual to leave \code{normalize.method="none"} in voom, i.e., the logCPM values should not usually be re-normalized in the \code{voom} call.

The \code{voom} method is similar in purpose to the limma-trend method, which uses \code{\link{eBayes}} or \code{\link{treat}} with \code{trend=TRUE}.
The voom method incorporates the mean-variance trend into the precision weights, whereas limma-trend incorporates the trend into the empirical Bayes moderation.
The voom method takes into account the sequencing depths (library sizes) of the individual columns of \code{counts} and applies the mean-variance trend on an individual observation basis.
limma-trend, on the other hand, assumes that the library sizes are not wildly different and applies the mean-variance trend on a genewise basis.
As noted by Law et al (2014), voom should be more powerful than limma-trend if the library sizes are very different but, otherwise, the two methods should give similar results.
}

\value{
An \code{\link[limma:EList]{EList}} object with the following components:
\item{E}{numeric matrix of normalized expression values on the log2 scale}
\item{weights}{numeric matrix of inverse variance weights}
\item{design}{design matrix}
\item{lib.size}{numeric vector of total normalized library sizes}
\item{genes}{dataframe of gene annotation extracted from \code{counts}}
\item{voom.xy}{if \code{save.plot}, list containing x and y coordinates for points in mean-variance plot}
\item{voom.line}{if \code{save.plot}, list containing coordinates of loess line in the mean-variance plot}
 }

\author{Charity Law and Gordon Smyth}

\references{
Law, CW, Chen, Y, Shi, W, Smyth, GK (2014).
Voom: precision weights unlock linear model analysis tools for RNA-seq read counts.
\emph{Genome Biology} 15, R29.
See also the Preprint Version at \url{http://www.statsci.org/smyth/pubs/VoomPreprint.pdf} incorporating some notational corrections.

Law, CW, Alhamdoosh, M, Su, S, Smyth, GK, Ritchie, ME (2016).
RNA-seq analysis is easy as 1-2-3 with limma, Glimma and edgeR.
\emph{F1000Research} 5, 1408.
\url{https://f1000research.com/articles/5-1408}

Law, CW, Alhamdoosh, M, Su, S, Dong, X, Tian, L, Smyth, GK, Ritchie, ME (2018).
RNA-seq analysis is easy as 1-2-3 with limma, Glimma and edgeR.
\emph{Bioconductor Workflow Package}.
\url{https://www.bioconductor.org/packages/RNAseq123/}
}

\seealso{
\code{\link{eBayes}},
\code{\link{voomWithQualityWeights}}.
\code{\link{vooma}} is similar to \code{voom} but for microarrays instead of RNA-seq.

A summary of functions for RNA-seq analysis is given in \link{11.RNAseq}.
}
\examples{
\dontrun{
keep <- filterByExpr(counts, design)
v <- voom(counts[keep,], design, plot=TRUE)
fit <- lmFit(v, design)
fit <- eBayes(fit, robust=TRUE)}
}
\keyword{rna-seq}
