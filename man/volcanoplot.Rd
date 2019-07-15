\title{Volcano Plot}
\name{volcanoplot}
\alias{volcanoplot}
\description{
Creates a volcano plot for a specified coefficient of a linear model.
}
\usage{
volcanoplot(fit, coef = 1, style = "p-value", highlight = 0, names = fit$genes$ID, hl.col="blue",
            xlab = "Log2 Fold Change", ylab = NULL, pch=16, cex=0.35, \dots)
}
\arguments{
  \item{fit}{an \code{MArrayLM} fitted linear model object.}
  \item{coef}{index indicating which coefficient of the linear model is to be plotted.}
  \item{style}{character string indicating which significance statistic to plot on the y-axis. Possibilities are \code{"p-value"} or \code{"B-statistic"}.}
  \item{highlight}{number of top genes to be highlighted by name.}
  \item{names}{character vector of length \code{nrow(fit)} giving gene names. Only used if \code{highlight > 0}.}
  \item{hl.col}{color for the gene names. Only used if \code{highlight > 0}.}
  \item{xlab}{character string giving label for x-axis}
  \item{ylab}{character string giving label for y-axis}
  \item{pch}{vector or list of plotting characters.}
  \item{cex}{numeric vector of plot symbol expansions.}
  \item{\dots}{any other arguments are passed to \code{plot}}
}

\details{
A volcano plot displays log fold changes on the x-axis versus a measure of statistical significance on the y-axis.
Here the significance measure can be -log(p-value) or the B-statistics, which give the posterior log-odds of differential expression.

The plot is optionally annotated with the names of the most significant genes.
}

\value{No value is returned but a plot is created on the current graphics device.}
\author{Gordon Smyth}
\seealso{
An overview of presentation plots following the fitting of a linear model in LIMMA is given in \link{06.LinearModels}.
}
\examples{
#  See lmFit examples
}
