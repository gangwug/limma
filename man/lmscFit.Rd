\name{lmscFit}
\alias{lmscFit}
\title{Fit Linear Model to Individual Channels of Two-Color Data}
\description{Fit a linear model to the individual log-intensities for each gene given a series of two-color arrays}
\usage{
lmscFit(object, design, correlation)
}
\arguments{
  \item{object}{an \code{\link[limma:malist]{MAList}} object or a list from which \code{M} and \code{A} values may be extracted}
  \item{design}{a numeric matrix containing the design matrix for linear model in terms of the individual channels.
  The number of rows should be twice the number of arrays.
  The number of columns will determine the number of coefficients estimated for each gene.}
  \item{correlation}{numeric value giving the intra-spot correlation}
}
\value{
An object of class \code{\link[limma:marraylm]{MArrayLM}}
}
\details{
For two color arrays, the channels measured on the same set of arrays are correlated.
The \code{M} and \code{A} however are uncorrelated for each gene.
This function fits a linear model to the set of M and A-values for each gene after re-scaling the M and A-values to have equal variances.
The input correlation determines the scaling required.
The input correlation is usually estimated using \code{\link{intraspotCorrelation}} before using \code{lmscFit}.

Missing values in \code{M} or \code{A} are not allowed.
}
\author{Gordon Smyth}
\seealso{
\code{\link[stats:lmfit]{lm.fit}}.

An overview of methods for single channel analysis in limma is given by \link{07.SingleChannel}.
}

\references{
Smyth, GK (2005). Individual channel analysis of two-colour microarray data.
\emph{Proceedings of the 55th Session of the International Statistics Institute}, 5-12 April 2005, Sydney, Australia;
Internatational Statistics Institute; Paper 116.
\url{http://www.statsci.org/smyth/pubs/ISI2005-116.pdf}

Smyth, GK, and Altman, NS (2013).
Separate-channel analysis of two-channel microarrays: recovering inter-spot information.
\emph{BMC Bioinformatics} 14, 165.
\url{http://www.biomedcentral.com/1471-2105/14/165}
}

\examples{
\dontrun{
# Subset of data from ApoAI case study in Limma User's Guide
# Avoid non-positive intensities
RG <- backgroundCorrect(RG,method="normexp")
MA <- normalizeWithinArrays(RG)
MA <- normalizeBetweenArrays(MA,method="Aq")
targets <- data.frame(Cy3=I(rep("Pool",6)),Cy5=I(c("WT","WT","WT","KO","KO","KO")))
targets.sc <- targetsA2C(targets)
targets.sc$Target <- factor(targets.sc$Target,levels=c("Pool","WT","KO"))
design <- model.matrix(~Target,data=targets.sc)
corfit <- intraspotCorrelation(MA,design)
fit <- lmscFit(MA,design,correlation=corfit$consensus)
cont.matrix <- cbind(KOvsWT=c(0,-1,1))
fit2 <- contrasts.fit(fit,cont.matrix)
fit2 <- eBayes(fit2)
topTable(fit2,adjust="fdr")
}
}
\keyword{models}
\keyword{regression}
