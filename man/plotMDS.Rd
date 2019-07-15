\title{Multidimensional scaling plot of distances between gene expression profiles}
\name{plotMDS}
\alias{plotMDS}
\alias{plotMDS.MDS}
\alias{plotMDS.default}
\alias{MDS-class}
\alias{show,MDS-method}

\description{
Plot samples on a two-dimensional scatterplot so that distances on the plot approximate the typical log2 fold changes between the samples.
}

\usage{
\method{plotMDS}{default}(x, top = 500, labels = NULL, pch = NULL, cex = 1,
     dim.plot = c(1,2), ndim = max(dim.plot), gene.selection = "pairwise",
     xlab = NULL, ylab = NULL, plot = TRUE, \dots)
\method{plotMDS}{MDS}(x, labels = NULL, pch = NULL, cex = 1, dim.plot = NULL,
     xlab = NULL, ylab = NULL, \dots)
}

\arguments{
  \item{x}{any data object which can be coerced to a matrix, for example an \code{ExpressionSet} or an \code{EList}.}
  \item{top}{number of top genes used to calculate pairwise distances.}
  \item{labels}{character vector of sample names or labels. Defaults to \code{colnames(x)}.}
  \item{pch}{plotting symbol or symbols. See \code{\link{points}} for possible values. Ignored if \code{labels} is non-\code{NULL}.}
  \item{cex}{numeric vector of plot symbol expansions.}
  \item{dim.plot}{integer vector of length two specifying which principal components should be plotted.}
  \item{ndim}{number of dimensions in which data is to be represented.}
  \item{gene.selection}{character, \code{"pairwise"} to choose the top genes separately for each pairwise comparison between the samples or \code{"common"} to select the same genes for all comparisons.}
  \item{xlab}{title for the x-axis.}
  \item{ylab}{title for the y-axis.}
  \item{plot}{logical. If \code{TRUE} then a plot is created on the current graphics device.}
  \item{\dots}{any other arguments are passed to \code{plot}, and also to \code{text} (if \code{pch} is \code{NULL}).}
}

\details{
This function is a variation on the usual multdimensional scaling (or principle coordinate) plot, in that a distance measure particularly appropriate for the microarray context is used.
The distance between each pair of samples (columns) is the root-mean-square deviation (Euclidean distance) for the top \code{top} genes.
Distances on the plot can be interpreted as \emph{leading log2-fold-change}, meaning
the typical (root-mean-square) log2-fold-change between the samples for the genes that distinguish those samples.

If \code{gene.selection} is \code{"common"}, then the top genes are those with the largest standard deviations between samples.
If \code{gene.selection} is \code{"pairwise"}, then a different set of top genes is selected for each pair of samples.
The pairwise feature selection may be appropriate for microarray data when different molecular pathways are relevant for distinguishing different pairs of samples.

If \code{pch=NULL}, then each sample is represented by a text label, defaulting to the column names of \code{x}.
If \code{pch} is not \code{NULL}, then plotting symbols are used.

See \code{\link[graphics]{text}} for possible values for \code{col} and \code{cex}.
}

\value{
If \code{plot=TRUE}, a plot is created on the current graphics device.

An object of class \code{"MDS"} is also invisibly returned.
This is a list containing the following components:
\item{distance.matrix}{numeric matrix of pairwise distances between columns of \code{x}}
\item{cmdscale.out}{output from the function \code{cmdscale} given the distance matrix}
\item{dim.plot}{dimensions plotted}
\item{x}{x-xordinates of plotted points}
\item{y}{y-cordinates of plotted points}
\item{gene.selection}{gene selection method}
}

\author{Di Wu and Gordon Smyth}

\references{
Ritchie, ME, Phipson, B, Wu, D, Hu, Y, Law, CW, Shi, W, and Smyth, GK (2015).
limma powers differential expression analyses for RNA-sequencing and microarray studies.
\emph{Nucleic Acids Research} 43, e47.
\url{http://nar.oxfordjournals.org/content/43/7/e47}
}

\seealso{
\code{\link{cmdscale}}

An overview of diagnostic functions available in LIMMA is given in \link{09.Diagnostics}.
}

\examples{
# Simulate gene expression data for 1000 probes and 6 microarrays.
# Samples are in two groups
# First 50 probes are differentially expressed in second group
sd <- 0.3*sqrt(4/rchisq(1000,df=4))
x <- matrix(rnorm(1000*6,sd=sd),1000,6)
rownames(x) <- paste("Gene",1:1000)
x[1:50,4:6] <- x[1:50,4:6] + 2
# without labels, indexes of samples are plotted.
mds <- plotMDS(x,  col=c(rep("black",3), rep("red",3)) )
# or labels can be provided, here group indicators:
plotMDS(mds,  col=c(rep("black",3), rep("red",3)), labels= c(rep("Grp1",3), rep("Grp2",3)))
}
