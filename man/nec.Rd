\name{nec}
\alias{nec}
\alias{neqc}
\title{NormExp Background Correction and Normalization Using Control Probes}
\description{Perform normexp background correction using negative control probes and quantile normalization using negative and positive control probes.
Particularly useful for Illumina BeadChips.}
\usage{
nec(x, status=NULL, negctrl="negative", regular="regular", offset=16,
    robust=FALSE, detection.p="Detection")
neqc(x, status=NULL, negctrl="negative", regular="regular", offset=16,
    robust=FALSE, detection.p="Detection", \dots)
}
\arguments{
  \item{x}{object of class \code{EListRaw} or \code{matrix} containing raw intensities for regular and control probes from a series of microarrays.}
  \item{status}{character vector giving probe types.  Defaults to \code{x$genes$Status} if \code{x} is an \code{EListRaw} object.}
  \item{negctrl}{character string identifier for negative control probes.}
  \item{regular}{character string identifier for regular probes, i.e., all probes other than control probes.}
  \item{offset}{numeric value added to the intensities after background correction.}
  \item{robust}{logical. Should robust estimators be used for the background mean and standard deviation?}
  \item{detection.p}{dection p-values.  Only used when no negative control probes can be found in the data.  Can be a numeric matrix or a character string giving the name of the component of \code{x$other} containing the matrix.}
  \item{\dots}{any other arguments are passed to \code{normalizeBetweenArrays.}}
  }
\details{
\code{neqc} performs background correction followed by quantile normalization, using negative control probes for background correction and both negative and positive controls for normalization (Shi et al, 2010).
\code{nec} is similar but performs background correction only.
These methods are particularly designed for Illumina BeadChip microarrays, but could be useful for other platforms for which good quality negative control probes or detection p-values are available.

When control data are available, these function call \code{\link{normexp.fit.control}} to estimate the parameters required by normal+exponential(normexp) convolution model with the help of negative control probes, followed by \code{\link{normexp.signal}} to perform the background correction.
If \code{x} contains background intensities \code{x$Eb}, then these are first subtracted from the foreground intensities, prior to normexp background correction.
After background correction, an \code{offset} is added to the data.

When expression values for negative controls are not available, the \code{detection.p} argument is used instead,
In that case, these functions call \code{\link{normexp.fit.detection.p}}, which infers the negative control probe intensities from the detection p-values associated with the regular probes.
The function outputs a message if this is done.

For more detailed descriptions of the arguments \code{x}, \code{status}, \code{negctrl}, \code{regular} and \code{detection.p}, please refer to functions \code{\link{normexp.fit.control}}, \code{\link{normexp.fit.detection.p}} and \code{\link{read.ilmn}}.

Both \code{nec} and \code{neqc} perform the above steps.
\code{neqc} continues on to quantile normalize the background-corrected intensities, including control probes.
After normalization, the intensities are log2 transformed and the control probes are removed.
}
\value{
\code{nec} produces a \code{\link{EListRaw-class}} or matrix object of the same dimensions as \code{x} containing background-corrected intensities, on the raw scale.
\code{neqc} produces a \code{\link{EList-class}} or matrix object containing normalized log2 intensities, with rows corresponding to control probes removed.
}

\references{
Shi W, Oshlack A and Smyth GK (2010).
Optimizing the noise versus bias trade-off for Illumina Whole Genome Expression BeadChips.
\emph{Nucleic Acids Research} 38, e204.
\url{http://nar.oxfordjournals.org/content/38/22/e204}
}

\author{Wei Shi and Gordon Smyth}

\seealso{ 
  An overview of background correction functions is given in \link{04.Background}.
  
  An overview of LIMMA functions for normalization is given in \link{05.Normalization}.

  \code{\link{normexp.fit.control}} estimates the parameters in the normal+exponential convolution model using the negative control probes.

  \code{\link{normexp.fit.detection.p}} estimates the parameters in the normal+exponential convolution model using negative control probe intensities inferred from regular probes by using their detection p values information.

  \code{\link{normexp.fit}} estimates parameters in the normal+exponential convolution model using a saddle-point approximation or other methods.

  \code{\link{neqc}} performs normexp background correction and quantile normalization aided by control probes.  
}

\examples{
\dontrun{
# neqc normalization for data which include control probes
x <- read.ilmn(files="sample probe profile.txt", ctrlfiles="control probe profile.txt")
y <- neqc(x)
fit <- lmFit(y,design)

# Same thing but in separate steps:
x.b <- nec(x)
y <- normalizeBetweenArrays(x.b,method="quantile")
y <- y[y$genes$Status=="regular",]

# neqc normalization for data without control probes
# neqc can process detection p-values in lieu of control probes
xr <- read.ilmn(files="sample probe profile.txt")
yr <- neqc(xr)
}
}

\keyword{background correction}
\keyword{illumina beadchips}
