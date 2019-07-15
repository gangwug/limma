\name{normalizeBetweenArrays}
\alias{normalizeBetweenArrays}
\title{Normalize Between Arrays}

\description{
Normalizes expression intensities so that the intensities or log-ratios have similar distributions across a set of arrays.
}

\usage{
normalizeBetweenArrays(object, method=NULL, targets=NULL, cyclic.method="fast", \dots)
}

\arguments{
  \item{object}{a numeric \code{matrix}, \code{\link[limma:EList]{EListRaw}}, \code{\link[limma:rglist]{RGList}} or \code{\link[limma:malist]{MAList}} object containing un-normalized expression data.
  If a matrix, then it is assumed to contain log-transformed single-channel data.}
  \item{method}{character string specifying the normalization method to be used.
  Choices for single-channel data are \code{"none"}, \code{"scale"}, \code{"quantile"} or \code{"cyclicloess"}.
  Choices for two-color data are those previously mentioned plus \code{"Aquantile"}, \code{"Gquantile"}, \code{"Rquantile"} or \code{"Tquantile"}.
  A partial string sufficient to uniquely identify the choice is permitted.
  The default is \code{"Aquantile"} for two-color data objects or \code{"quantile"} for single-channel objects.}
  \item{targets}{vector, factor or matrix of length twice the number of arrays, used to indicate target groups if \code{method="Tquantile"}}
  \item{cyclic.method}{character string indicating the variant of \code{normalizeCyclicLoess} to be used if \code{method=="cyclicloess"}, see \code{\link{normalizeCyclicLoess}} for possible values.}
  \item{\dots}{other arguments are passed to \code{normalizeQuantiles} or \code{normalizeCyclicLoess}}
}

\details{
\code{normalizeBetweenArrays} normalizes expression values to achieve consistency between arrays.
For two-color arrays, normalization between arrays is usually a follow-up step after normalization within arrays using \code{normalizeWithinArrays}.
For single-channel arrays, within array normalization is not usually relevant and so \code{normalizeBetweenArrays} is the sole normalization step.

For single-channel data, the scale, quantile or cyclic loess normalization methods can be applied to the columns of data.
Trying to apply other normalization methods when \code{object} is a \code{matrix} or \code{EListRaw} object will produce an error.
If \code{object} is an \code{EListRaw} object, then normalization will be applied to the matrix \code{object$E} of expression values, which will then be log2-transformed.
Scale (\code{method="scale"}) scales the columns to have the same median.
Quantile and cyclic loess normalization was originally proposed by Bolstad et al (2003) for Affymetrix-style single-channel arrays.
Quantile normalization forces the entire empirical distribution of each column to be identical.
Cyclic loess normalization applies loess normalization to all possible pairs of arrays, usually cycling through all pairs several times.
Cyclic loess is slower than quantile, but allows probe-wise weights and is more robust to unbalanced differential expression.

The other normalization methods are for two-color arrays.
Scale normalization was proposed by Yang et al (2001, 2002) and is further explained by Smyth and Speed (2003).
The idea is simply to scale the log-ratios to have the same median-absolute-deviation (MAD) across arrays.
This idea has also been implemented by the \code{maNormScale} function in the marray package.
The implementation here is slightly different in that the MAD scale estimator is replaced with the median-absolute-value and the A-values are normalized as well as the M-values.

Quantile normalization was explored by Yang and Thorne (2003) for two-color cDNA arrays.
\code{method="quantile"} ensures that the intensities have the same empirical distribution across arrays and across channels.
\code{method="Aquantile"} ensures that the A-values (average intensities) have the same empirical distribution across arrays leaving the M-values (log-ratios) unchanged.
These two methods are called "q" and "Aq" respectively in Yang and Thorne (2003).

\code{method="Tquantile"} performs quantile normalization separately for the groups indicated by \code{targets}.
\code{targets} may be a target frame such as read by \code{readTargets} or can be a vector indicating green channel groups followed by red channel groups.

\code{method="Gquantile"} ensures that the green (first) channel has the same empirical distribution across arrays, leaving the M-values (log-ratios) unchanged.
This method might be used when the green channel is a common reference throughout the experiment.
In such a case the green channel represents the same target throughout, so it makes compelling sense to force the distribution of intensities to be same for the green channel on all the arrays, and to adjust to the red channel accordingly.
\code{method="Rquantile"} ensures that the red (second) channel has the same empirical distribution across arrays, leaving the M-values (log-ratios) unchanged.
Both \code{Gquantile} and \code{Rquantile} normalization have the implicit effect of changing the red and green log-intensities by equal amounts.

See the limma User's Guide for more examples of use of this function.
}

\value{
If \code{object} is a matrix then \code{normalizeBetweenArrays} produces a matrix of the same size.
If \code{object} is an \code{EListRaw} object, then an \code{EList} object with expression values on the log2 scale is produced.
For two-color data, \code{normalizeBetweenArrays} produces an \code{\link[limma:malist]{MAList}} object with M and A-values on the log2 scale.
}

\author{Gordon Smyth}

\references{
Bolstad, B. M., Irizarry R. A., Astrand, M., and Speed, T. P. (2003), A comparison of normalization methods for high density oligonucleotide array data based on bias and variance. \emph{Bioinformatics} \bold{19}, 185-193.

Smyth, G. K., and Speed, T. P. (2003). Normalization of cDNA microarray data. \emph{Methods} \bold{31}, 265-273. 

Yang, Y. H., Dudoit, S., Luu, P., and Speed, T. P. (2001). Normalization for cDNA microarray data. In \emph{Microarrays: Optical Technologies and Informatics}, M. L. Bittner, Y. Chen, A. N. Dorsel, and E. R. Dougherty (eds), Proceedings of SPIE, Volume 4266, pp. 141-152. 

Yang, Y. H., Dudoit, S., Luu, P., Lin, D. M., Peng, V., Ngai, J., and Speed, T. P. (2002). Normalization for cDNA microarray data: a robust composite method addressing single and multiple slide systematic variation. \emph{Nucleic Acids Research} \bold{30}(4):e15.

Yang, Y. H., and Thorne, N. P. (2003). Normalization for two-color cDNA microarray data.
In: D. R. Goldstein (ed.), \emph{Science and Statistics: A Festschrift for Terry Speed}, IMS Lecture Notes - Monograph Series, Volume 40, pp. 403-418.
}

\seealso{
  An overview of LIMMA functions for normalization is given in \link{05.Normalization}.

  The \code{\link{neqc}} function provides a variation of quantile normalization that is customized for Illumina BeadChips.
  This method uses control probes to refine the background correction and normalization steps.

  Note that vsn normalization, previously offered as a method of this function, is now performed by the \code{\link{normalizeVSN}} function.

  See also \code{maNormScale} in the marray package and
  \code{\link[affy:normalize-methods]{normalize-methods}} in the affy package.
}

\examples{
ngenes <- 100
narrays <- 4
x <- matrix(rnorm(ngenes*narrays),100,4)
y <- normalizeBetweenArrays(x)
}

\keyword{normalization}
