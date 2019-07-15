\name{05.Normalization}
\alias{05.Normalization}
\title{Topic: Normalization of Microarray Data}

\description{
This page gives an overview of the LIMMA functions available to normalize data from single-channel or two-colour microarrays.
Smyth and Speed (2003) give an overview of the normalization techniques implemented in the functions for two-colour arrays.

Usually data from spotted microarrays will be normalized using \code{\link{normalizeWithinArrays}}.
A minority of data will also be normalized using \code{\link{normalizeBetweenArrays}} if diagnostic plots suggest a difference in scale between the arrays.

In rare circumstances, data might be normalized using \code{\link{normalizeForPrintorder}} before using \code{\link{normalizeWithinArrays}}.

All the normalization routines take account of spot quality weights which might be set in the data objects.
The weights can be temporarily modified using \code{\link{modifyWeights}} to, for example, remove ratio control spots from the normalization process.

If one is planning analysis of single-channel information from the microarrays rather than analysis of differential expression based on log-ratios, then the data should be normalized using a single channel-normalization technique.
Single channel normalization uses further options of the \code{\link{normalizeBetweenArrays}} function.
For more details see the \link[=limmaUsersGuide]{LIMMA User's Guide} which includes a section on single-channel normalization.

\code{normalizeWithinArrays} uses utility functions \code{\link{MA.RG}}, \code{\link{loessFit}} and \code{\link{normalizeRobustSpline}}.

\code{normalizeBetweenArrays} is the main normalization function for one-channel arrays,
as well as an optional function for two-colour arrays.
\code{normalizeBetweenArrays} uses utility functions \code{\link{normalizeMedianAbsValues}}, \code{\link{normalizeMedianAbsValues}}, \code{\link{normalizeQuantiles}} and \code{\link{normalizeCyclicLoess}}, none of which need to be called directly by users.

\code{\link{neqc}} is a between array normalization function customized for Illumina BeadChips.

The function \code{\link{normalizeVSN}} is also provided as a interface to the vsn package.
It performs variance stabilizing normalization, an algorithm which includes background correction, within and between normalization together, and therefore doesn't fit into the paradigm of the other methods.

\code{\link{removeBatchEffect}} can be used to remove a batch effect, associated with hybridization time or some other technical variable, prior to unsupervised analysis.
}

\author{Gordon Smyth}
\references{
Smyth, G. K., and Speed, T. P. (2003). Normalization of cDNA microarray data.
\emph{Methods} 31, 265-273.
\url{http://www.statsci.org/smyth/pubs/normalize.pdf}
}

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
