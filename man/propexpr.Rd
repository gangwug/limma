\name{propexpr}
\alias{propexpr}
\title{Estimate Proportion of Expressed Probes}
\description{Estimate the proportion of microarray probes which are expressed in each array.}
\usage{
propexpr(x, neg.x=NULL, status=x$genes$Status, labels=c("negative","regular"))
}
\arguments{
  \item{x}{matrix or similar object containing raw intensities for a set of arrays.}
  \item{neg.x}{matrix or similar object containing raw intensities for negative control probes for the same arrays. If \code{NULL}, then negative controls must be provided in \code{x}.}
  \item{status}{character vector specifying control type of each probe. Only used if \code{neg.x} is \code{NULL}.}
  \item{labels}{character vector giving the \code{status} values for negative control probes and regular (non-control) probes respectively. If of length 1, then all probes other than the negative controls are assumed to be regular. Only used if \code{neg.x} is \code{NULL}.}
  }

\details{
This function estimates the overall proportion of probes on each microarray that are correspond to expressed genes using the method of Shi et al (2010).
The function is especially useful for Illumina BeadChips arrays, although it can in principle be applied to any platform with good quality negative controls.

The negative controls can be supplied either as rows of \code{x} or as a separate matrix.
If supplied as rows of \code{x}, then the negative controls are identified by the \code{status} vector.
\code{x} might also include other types of control probes, but these will be ignored in the calculation.

Illumina BeadChip arrays contain 750~1600 negative control probes.
If \code{read.idat} is used to read Illumina expression IDAT files, then the control probes will be populated as rows of the output \code{EListRaw} object, and the vector \code{x$genes$Status} will be set to identify control probes.

Alternatively, expression values can be exported from Illumina's GenomeStudio software as tab-delimited text files.
In this case, the control probes are usually written to a separate file from the regular probes.
}

\value{
Numeric vector giving the proportions of expressed probes in each array.
}

\references{
Shi, W, de Graaf, C, Kinkel, S, Achtman, A, Baldwin, T, Schofield, L, Scott, H, Hilton, D, Smyth, GK (2010).
Estimating the proportion of microarray probes expressed in an RNA sample.
\emph{Nucleic Acids Research} 38(7), 2168-2176.
\url{https://www.ncbi.nlm.nih.gov/pubmed/20056656}
}

\author{Wei Shi and Gordon Smyth}

\seealso{ 
Description to the control probes in Illumina BeadChips can be found in \code{\link{read.ilmn}}.
}

\examples{
\dontrun{
# Read Illumina binary IDAT files
x <- read.idat(idat, bgx)
propexpr(x)

# Read text files exported from GenomeStudio
x <- read.ilmn(files = "sample probe profile.txt",
               ctrlfiles = "control probe profile.txt")
propexpr(x)
}
}

\keyword{background correction}
\keyword{illumina beadchips}
