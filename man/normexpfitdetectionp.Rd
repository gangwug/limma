\name{normexp.fit.detection.p}
\alias{normexp.fit.detection.p}
\title{Estimate Normexp Model Parameter Using Negative Controls Inferred from Regular Probes}
\description{Detection p values from Illumina BeadChip microarray data can be used to infer negative control probe intensities from regular probe intensities by using detection p value information when negative control data are not available. The inferred negative control intensities can then be used in the background correction in the same way as those control data outputted from BeadChip used in the \code{\link{normexp.fit.control}} function.
}
\usage{
normexp.fit.detection.p(x, detection.p="Detection")
}
\arguments{
  \item{x}{object of class \code{EListRaw-class} or \code{matrix} containing raw intensities of regular probes for a series of microarrays}
  \item{detection.p}{a character string giving the name of the component which contains detection p value information in \code{x} or a numeric matrix giving detection p values, \code{Detection} by default}
  }
\details{
This function estimates the normexp parameters in the same way as \code{\link{normexp.fit.control}} does, except that negative control probe intensities are inferred from regular probes by taking advantage of detection p value information rather than from the control probe profile outputted by BeadStudio.

Calculation of detection p values in Illumina BeadChip data is based on the rank of probe intensities in the list of negative control probe intensities.
Therefore, the detection p values can be used to find regular probes which have expression intensities falling into the range of negative control probe intensities. 
These probes give a good approximation to the real negative control data and thus can be used to estimate the mean and standard deviation of background intensities when negative control data is not available.

If \code{x} is an \code{EListRaw-class} object, this function will try to look for the component which includes detection p value matrix in \code{x} when \code{detection.p} is a character string. 
This function assumes that this component is located within the \code{other} component in \code{x}.
The component name specified by \code{detection.p} should be exactly the same as the name of the detection p value component in \code{x}.
If \code{detection.p} is a matrix, then this matrix will be used as the detection p value data used in this function.

If \code{x} is an \code{matrix} object, then \code{detection.p} has to be a data matrix which includes detection p values.

When \code{detection.p} is a \code{matrix}, it has to have the same dimension as that of \code{x}.

This function will replace the detection p values with 1 subtracted by these values if high intensity probes have detection p values less than those from low intensity probes.

Note that when control data are available, the \code{\link{normexp.fit.control}} function should be used instead.
}
\value{
A matrix containing estimated parameters with rows being arrays and with columns being parameters.
Column names are \code{mu}, \code{logsigma} and \code{logalpha}.
}

\references{
Shi W, Oshlack A and Smyth GK (2010).
Optimizing the noise versus bias trade-off for Illumina Whole Genome Expression BeadChips.
\emph{Nucleic Acids Research} 38, e204.
\url{http://nar.oxfordjournals.org/content/38/22/e204}
}

\author{Wei Shi and Gordon Smyth}

\seealso{
\code{\link{nec}} calls this function to get the parameters of the normal+exponential convolution model when control probe profile file is not available and then calls \code{\link{normexp.signal}} to perform the background correction.

\code{\link{normexp.fit.control}} estimates normexp parameters using control data outputted by BeadStudio.

\code{\link{normexp.fit}} estimates normexp parameters using a saddle-point approximation or other mothods.

An overview of background correction functions is given in \code{\link{04.Background}}.
}

\examples{
\dontrun{
# read in BeadChip data which do not have control data available
x <- read.ilmn(files="sample probe profile")
# estimated normexp parameters
normexp.fit.detection.p(x)
# normalization using inferred negative controls
y <- neqc(x)
}
}

\keyword{models}
