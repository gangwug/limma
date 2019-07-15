\name{beadCountWeights}
\alias{beadCountWeights}

\title{Bead Count Weights for Illumina BeadChips}

\description{
Estimates weights which account for biological variation and technical variation resulting from varying bead numbers.
}

\usage{
beadCountWeights(y, x, design = NULL, bead.stdev = NULL, bead.stderr = NULL,
                 nbeads = NULL, array.cv = TRUE, scale = FALSE)
}

\arguments{
 \item{y}{an \code{"EList"} object or a numeric matrix containing normalized log2-expression values.}
 \item{x}{an \code{"EListRaw"} object or a numeric matrix of raw expression values, with same dimensions as \code{y}.}
 \item{design}{the design matrix of the microarray experiment, with rows
           corresponding to arrays and columns to coefficients to be
           estimated.  Defaults to \code{y$design} or, if that is \code{NULL}, then to a column of ones meaning that the
           arrays are treated as replicates.}
 \item{bead.stdev}{numeric matrix of bead-level standard deviations.}
 \item{bead.stderr}{numeric matrix of bead-level standard errors. Not required if \code{bead.stdev} is set.}
 \item{nbeads}{numeric matrix containing number of beads.}
 \item{array.cv}{logical, should technical variation for each observation be calculated from a constant or array-specific coefficient of variation?  The default is to use array-specific coefficients of variation.}
 \item{scale}{logical, should weights be scaled so that the average weight size is the mean of the inverse technical variance along a probe? By default, weights are scaled so that the average weight size along a probe is 1.}
}

\details{
This function estimates optimum weights using the bead statistics for each probe for an Illumina expression BeadChip.
It can be used with any Illumina expression BeadChip, but is most likely to be useful with HumanHT-12 BeadChips.

Arguments \code{x} and \code{y} are both required.
\code{x} contains the raw expression values and \code{y} contains the corresponding log2 values for the same probes and the same arrays after background correction and normalization.
\code{x} and \code{y} be any type of object that can be coerced to a matrix, with rows corresponding to probes and columns to arrays.
\code{x} and \code{y} must contain the same rows and columns in the same order.

The reliability of the normalized expression value for each probe on each array is measured by estimating its technical and biological variability.
The bead number weights are the inverse sum of the technical and biological variances.

The technical variance for each probe on each array is inversely proportional to the number of beads and is estimated using array-specific bead-level coefficients of variation.

Coefficients of variation are calculated using raw expression values. 

The biological variance for each probe across the arrays are estimated using a Newton iteration, with the assumption that the total residual deviance for each probe from \code{lmFit} is inversely proportional to the sum of the technical variance and biological variance. 

Only one of \code{bead.stdev} or \code{bead.stderr} needs to be set.
If \code{bead.stdev} is not provided, then it will be computed as \code{bead.stderr * sqrt(nbeads)}.

If arguments \code{bead.stdev} and \code{nbeads} are not set explicitly in the call, then they will be extracted from\code{y$other$BEAD_STDEV} and \code{y$other$Avg_NBEADS}.
An \code{EList} object containing these components can be created by \code{read.idat} or \code{read.ilmn}, see the example code below.
}

\value{
A list object with the following components:

\item{weights}{numeric matrix of bead number weights}
\item{cv.constant}{numeric value of constant bead-level coefficient of variation}
\item{cv.array}{numeric vector of array-specific bead-level coefficient of variation}
\item{var.technical}{numeric matrix of technical variances}
\item{var.biological}{numeric vector of biological variances}
}

\references{
Law, CW (2013).
\emph{Precision weights for gene expression analysis}.
PhD Thesis. University of Melbourne, Australia.
\url{http://repository.unimelb.edu.au/10187/17598}
}

\author{Charity Law and Gordon Smyth}

\examples{
\dontrun{
z <- read.ilmn(files="probesummaryprofile.txt",
     ctrfiles="controlprobesummary.txt",
     other.columns=c("BEAD_STDEV","Avg_NBEADS"))
y <- neqc(z)
x <- z[z$genes$Status=="regular",]
bcw <- beadCountWeights(y,x,design)
fit <- lmFit(y,design,weights=bcw$weights)
fit <- eBayes(fit)
}
}

\seealso{
\code{\link{read.ilmn}}, \code{\link{read.idat}}, \code{\link{neqc}}.

An overview of linear model functions in limma is given by \link{06.LinearModels}.
}
