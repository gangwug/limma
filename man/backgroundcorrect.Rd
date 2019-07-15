\name{backgroundCorrect}
\alias{backgroundCorrect}
\alias{backgroundCorrect.matrix}
\title{Correct Intensities for Background}
\description{
Background correct microarray expression intensities.
}
\usage{
backgroundCorrect(RG, method="auto", offset=0, printer=RG$printer,
                  normexp.method="saddle", verbose=TRUE)
backgroundCorrect.matrix(E, Eb=NULL, method="auto", offset=0, printer=NULL,
                  normexp.method="saddle", verbose=TRUE)
}
\arguments{
  \item{RG}{a numeric matrix, \code{\link[limma:EList]{EListRaw}} or \code{\link[limma:rglist]{RGList}} object.}
  \item{E}{numeric matrix containing foreground intensities.}
  \item{Eb}{numeric matrix containing background intensities.}
  \item{method}{character string specifying correction method.  Possible values are \code{"auto"}, \code{"none"}, \code{"subtract"}, \code{"half"}, \code{"minimum"}, \code{"movingmin"}, \code{"edwards"} or \code{"normexp"}.
  If \code{RG} is a matrix, possible values are restricted to \code{"none"} or \code{"normexp"}.
  The default \code{"auto"} is interpreted as \code{"subtract"} if background intensities are available or \code{"normexp"} if they are not.}
  \item{offset}{numeric value to add to intensities}
  \item{printer}{a list containing printer layout information, see \code{\link{PrintLayout-class}}. Ignored if \code{RG} is a matrix.}
  \item{normexp.method}{character string specifying parameter estimation strategy used by normexp, ignored for other methods. Possible values are \code{"saddle"}, \code{"mle"}, \code{"rma"} or \code{"rma75"}.}
  \item{verbose}{logical. If \code{TRUE}, progress messages are sent to standard output}
}
\details{
This function implements the background correction methods reviewed or developed in Ritchie et al (2007) and Silver at al (2009).
Ritchie et al (2007) recommend \code{method="normexp"} whenever \code{RG} contains local background estimates.
Silver et al (2009) shows that either \code{normexp.method="mle"} or \code{normexp.method="saddle"} are excellent options for normexp.
If \code{RG} contains morphological background estimates instead (available from SPOT or GenePix image analysis software), then \code{method="subtract"} performs well.

If \code{method="none"} then no correction is done, i.e., the background intensities are treated as zero.
If \code{method="subtract"} then the background intensities are subtracted from the foreground intensities.
This is the traditional background correction method, but is not necessarily recommended.
If \code{method="movingmin"} then the background estimates are replaced with the minimums of the backgrounds of the spot and its eight neighbors, i.e., the background is replaced by a moving minimum of 3x3 grids of spots.

The remaining methods are all designed to produce positive corrected intensities.
If \code{method="half"} then any intensity which is less than 0.5 after background subtraction is reset to be equal to 0.5.
If \code{method="minimum"} then any intensity which is zero or negative after background subtraction is set equal to half the minimum of the positive corrected intensities for that array.
If \code{method="edwards"} a log-linear interpolation method is used to adjust lower intensities as in Edwards (2003).
If \code{method="normexp"} a convolution of normal and exponential distributions is fitted to the foreground intensities using the background intensities as a covariate, and the expected signal given the observed foreground becomes the corrected intensity.
This results in a smooth monotonic transformation of the background subtracted intensities such that all the corrected intensities are positive.

The normexp method is available in a number of variants depending on how the model parameters are estimated, and these are selected by \code{normexp.method}.
Here \code{"saddle"} gives the saddle-point approximation to maximum likelihood from Ritchie et al (2007) and improved by Silver et al (2009), \code{"mle"} gives exact maximum likelihood from Silver at al (2009), \code{"rma"} gives the background correction algorithm from the RMA-algorithm for Affymetrix microarray data as implemented in the affy package, and \code{"rma75"} gives the RMA-75 method from McGee and Chen (2006).
In practice \code{"mle"} performs well and is nearly as fast as \code{"saddle"}, but \code{"saddle"} is the default for backward compatibility.
See \code{\link{normexp.fit}} for more details.

The \code{offset} can be used to add a constant to the intensities before log-transforming, so that the log-ratios are shrunk towards zero at the lower intensities.
This may eliminate or reverse the usual 'fanning' of log-ratios at low intensities associated with local background subtraction.

Background correction (background subtraction) is also performed by the \code{\link{normalizeWithinArrays}} method for \code{RGList} objects, so it is not necessary to call \code{backgroundCorrect} directly unless one wants to use a method other than simple subtraction.
Calling \code{backgroundCorrect} before \code{normalizeWithinArrays} will over-ride the default background correction.
}

\value{
A matrix, \code{EListRaw} or \code{RGList} object in which foreground intensities have been background corrected and any components containing background intensities have been removed.
}
\references{
Edwards, D. E. (2003). Non-linear normalization and background correction in one-channel cDNA microarray studies
\emph{Bioinformatics} 19, 825-833. 

McGee, M., and Chen, Z. (2006).
Parameter estimation for the exponential-normal convolution model for background correction of Affymetrix GeneChip data.
\emph{Stat Appl Genet Mol Biol}, Volume 5, Article 24.

Ritchie, M. E., Silver, J., Oshlack, A., Silver, J., Holmes, M., Diyagama, D., Holloway, A., and Smyth, G. K. (2007).
A comparison of background correction methods for two-colour microarrays.
\emph{Bioinformatics} 23, 2700-2707.
\url{http://bioinformatics.oxfordjournals.org/content/23/20/2700}

Silver, J., Ritchie, M. E., and Smyth, G. K. (2009).
Microarray background correction: maximum likelihood estimation for the normal-exponential convolution model.
\emph{Biostatistics} 10, 352-363.
\url{http://biostatistics.oxfordjournals.org/content/10/2/352}
}
\author{Gordon Smyth}
\examples{
RG <- new("RGList", list(R=c(1,2,3,4),G=c(1,2,3,4),Rb=c(2,2,2,2),Gb=c(2,2,2,2)))
backgroundCorrect(RG)
backgroundCorrect(RG, method="half")
backgroundCorrect(RG, method="minimum")
backgroundCorrect(RG, offset=5)
}
\seealso{
\code{\link{kooperberg}}, \code{\link{neqc}}.

An overview of background correction functions is given in \code{\link{04.Background}}.
}

\keyword{background correction}
