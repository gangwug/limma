\name{normexp.fit.control}
\alias{normexp.fit.control}
\title{Normexp Model Parameter Estimation Aided by Negative Controls}
\description{The mean and log-standard-deviation of the background-normal part of the normexp+exponential convolution model is estimated as the mean and log-standard deviation of intensities from negative control probes. The log-mean of the signal-exponential part is estimated as the log of the difference between signal mean and background mean.}
\usage{
normexp.fit.control(x, status=NULL, negctrl="negative", regular="regular", robust=FALSE)
}
\arguments{
  \item{x}{object of class \code{EListRaw-class} or \code{matrix} containing raw intensities for regular and control probes for a series of microarrays}
  \item{status}{character vector giving probe types.}
  \item{negctrl}{character string identifier for negative control probes.}
  \item{regular}{character string identifier for regular probes.}
  \item{robust}{logical. Should robust estimators be used for the background mean and standard deviation?}
  }
\details{
\code{x} has to contain raw expression intensities from both regular probes and negative control probes.

The probe type information for an object of \code{\link{EListRaw-class}} is normally saved in the \code{Status} column of its \code{genes} component.
However, it will be overriden by the \code{status} parameter if it is explicitly provided to this function.
If \code{x} is a \code{matrix} object, the probe type information has to be provided through the \code{status} parameter of this function.
Regular probes have the status \code{regular}.
Negative control probes have the status indicated by \code{negctrl}, which is \code{negative} by default.

This function estimates parameters of the normal+exponential convolution model with the help of negative control probes.
The mean and log-standard-deviation of the background-normal part of the normexp+exponential(normexp) convolution model are estimated as the mean and log-standard deviation of intensities from negative control probes respectively.
The log-mean of the signal-exponential part is estimated as the log of the difference between signal mean and background mean.
The signal mean is simply the mean of intensities from regular probes.

When negative control probes are not available, the \code{\link{normexp.fit.detection.p}} function can be used to estimate the normexp model parameters which infers the negative control probe intensities from regular probes by taking advantage of their detection p value information.
}
\value{
A matrix containing estimated parameters with rows being arrays and with columns being parameters.
Column names are \code{mu}, \code{logsigma} and \code{logalpha}.
}

\references{
Shi W, Oshlack A and Smyth GK (2010). Optimizing the noise versus bias trade-off for Illumina Whole Genome Expression BeadChips. Nucleic Acids Research, 38(22):e204. Epub 2010 Oct 6. PMID: 20929874 
}

\author{Wei Shi and Gordon Smyth}

\seealso{
\code{\link{nec}} calls this function to get the parameters of the normal+exponential convolution model and then calls \code{\link{normexp.signal}} to perform the background correction.

\code{\link{normexp.fit.detection.p}} estimates the parameters in the normal+exponential convolution model using negative control probe intensities inferred from regular probes by using their detection p values information.

\code{\link{normexp.fit}} estimates normexp parameters using a saddle-point approximation or other mothods.

An overview of background correction functions is given in \code{\link{04.Background}}.
}

\examples{
\dontrun{
# read in BeadChip probe profile file and control profile file
x <- read.ilmn(files="sample probe profile", ctrlfiles="control probe profile")
# estimated normexp parameters
normexp.fit.control(x)
# normalization using control data
y <- neqc(x)
}
}

\keyword{models}
