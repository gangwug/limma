\name{detectionPValues}
\alias{detectionPValues.default}
\alias{detectionPValues.EListRaw}
\alias{detectionPValues}

\title{Detection P-Values from Negative Controls}

\description{Compute the proportion of negative controls greater than each observed expression value.
Particularly useful for Illumina BeadChips.}

\usage{
\method{detectionPValues}{EListRaw}(x, status = NULL, \dots)
\method{detectionPValues}{default}(x, status, negctrl = "negative", \dots)
}

\arguments{
  \item{x}{object of class \code{EListRaw} or a numeric \code{matrix} containing raw intensities for regular and control probes from a series of microarrays.}
  \item{status}{character vector giving probe types.  Defaults to \code{x$genes$Status} if \code{x} is an \code{EListRaw} object.}
  \item{negctrl}{character string identifier for negative control probes.}
  \item{\dots}{other arguments are not currently used.}
  }

\details{
The rows of \code{x} for which \code{status == negctrl} are assumed to correspond to negative control probes.

For each column of \code{x}, the detection p-values are defined as \code{(N.eq/2 + N.gt) / N.neg}, where \code{N.gt} is the number of negative controls with expression greater than the observed value, \code{N.eq} is the number of negative controls with expression equal to the observed value, and \code{N.neg} is the total number of negative controls.

When used on Illumina BeadChip data, this function produces essentially the same detection p-values as returned by Illumina's GenomeStudio software.
}

\value{numeric matrix of same dimensions as \code{x} containing detection p-values.}

\references{
Shi, W, de Graaf, C, Kinkel, S, Achtman, A, Baldwin, T, Schofield, L, Scott, H, Hilton, D, Smyth, GK (2010).
Estimating the proportion of microarray probes expressed in an RNA sample.
\emph{Nucleic Acids Research} 38(7), 2168-2176.
\url{https://www.ncbi.nlm.nih.gov/pubmed/20056656}
}

\author{Gordon Smyth}

\seealso{ 
  An overview of LIMMA functions to read expression data is given in \link{03.ReadingData}.

  \code{\link{read.idat}} reads Illumina BeadChip expression data from binary IDAT files.

  \code{\link{neqc}} performs normexp background correction and quantile normalization aided by control probes.  
}

\examples{
\dontrun{
# Read Illumina binary IDAT files
x <- read.idat(idat, bgx)
x$genes$DectionPValue <- detectionPValues(x)
y <- neqc(x)
}
}

\keyword{background correction}
\keyword{illumina beadchips}
