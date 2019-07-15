\name{normalizeCyclicLoess}
\alias{normalizeCyclicLoess}
\title{Normalize Columns of a Matrix by Cyclic Loess}
\description{
Normalize the columns of a matrix, cyclicly applying loess normalization to normalize each pair of columns to each other.
}
\usage{
normalizeCyclicLoess(x, weights = NULL, span=0.7, iterations = 3, method = "fast")
}
\arguments{
  \item{x}{numeric matrix, or object which can be coerced to a numeric matrix, containing log-expression values.}
  \item{weights}{numeric vector of probe weights. Must be non-negative.}
  \item{span}{span of loess smoothing window, between 0 and 1.}
  \item{iterations}{number of times to cycle through all pairs of columns.}
  \item{method}{character string specifying which variant of the cyclic loess method to use. Options are \code{"fast"}, \code{"affy"} or \code{"pairs"}.}
}
\details{
This function is intended to normalize single channel or A-value microarray intensities between arrays.
Cyclic loess normalization is similar effect and intention to quantile normalization, but with some advantages, in particular the ability to incorporate probe weights.

A number of variants of cylic loess have been suggested.
\code{method="pairs"} implements the intuitive idea that each pair of arrays is subjected to loess normalization as for two-color arrays.
This process is simply cycled through all possible pairs of arrays, then repeated for several \code{iterations}.
This is the method described by Ballman et al (2004) as ordinary cyclic loess normalization.

\code{method="affy"} implements a method similar to \code{normalize.loess} in the affy package,
except that here we call \code{lowess} instead of \code{loess} and avoid the use of probe subsets and the \code{predict} function.
In this approach, no array is modified until a complete cycle of all pairs has been completed.
The adjustments are stored for a complete iteration, then averaged, and finally used to modify the arrays.
The \code{"affy"} method is invariant to the order of the columns of \code{x}, whereas the \code{"pairs"} method is not.
The affy approach is presumably that used by Bolstad et al (2003), although the algorithm was not explicitly described in that article.

\code{method="fast"} implements the "fast linear loess" method of Ballman et al (2004), whereby each array is simply normalized to a reference array,
the reference array being the average of all the arrays.
This method is relatively fast because computational time is linear in the number of arrays, whereas \code{"pairs"} and \code{"affy"} are quadratic in the number of arrays.
\code{"fast"} requires n lowess fits per iteration, where n is the number of arrays, whereas \code{"pairs"} and \code{"affy"} require n*(n-1)/2 lowess fits per iteration.
}
\value{
A matrix of the same dimensions as \code{x} containing the normalized values.
}
\references{
Bolstad, B. M., Irizarry R. A., Astrand, M., and Speed, T. P. (2003). A comparison of normalization methods for high density oligonucleotide array data based on bias and variance. \emph{Bioinformatics} \bold{19}, 185-193.

Ballman, KV
Grill, DE,
Oberg, AL and
Therneau, TM (2004).
Faster cyclic loess: normalizing RNA arrays via linear models.
\emph{Bioinformatics} 20, 2778-2786. 
}
\author{Yunshun (Andy) Chen and Gordon Smyth}
\seealso{
An overview of LIMMA functions for normalization is given in \link{05.Normalization}.

\link[affy]{normalize.loess} in the affy package also implements cyclic loess normalization, without weights.
}  

\keyword{normalization}
