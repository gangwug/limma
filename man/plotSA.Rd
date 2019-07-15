\title{Sigma vs A plot for microarray linear model}
\name{plotSA}
\alias{plotSA}
\description{
Plot residual standard deviation versus average log expression for a fitted microarray linear model.
}
\usage{
plotSA(fit, xlab = "Average log-expression", ylab = "sqrt(sigma)", zero.weights = FALSE,
       pch = 16, cex = 0.3, col = c("black","red"), \dots)
}
\arguments{
  \item{fit}{an \code{MArrayLM} object.}
  \item{xlab}{label for x-axis}
  \item{ylab}{label for y-axis}
  \item{zero.weights}{logical, should genes with all zero weights be plotted?}
  \item{pch}{vector of codes for plotting characters.}
  \item{cex}{numeric, vector of expansion factors for plotting characters.}
  \item{col}{plotting colors for regular and outlier variances respectively.}
  \item{\dots}{any other arguments are passed to \code{plot}}
}

\details{
This plot is used to check the mean-variance relationship of the expression data, after fitting a linear model.
A scatterplot of residual-variances vs average log-expression is created.
The plot is especially useful for examining the mean-variance trend estimated by \code{\link{eBayes}} or \code{\link{treat}} with \code{trend=TRUE}.
It can be considered as a routine diagnostic plot in the limma-trend pipeline.

If robust empirical Bayes was used to create \code{fit}, then outlier variances are highlighted in the color given by \code{col[2]}.

The y-axis is square-root \code{fit$sigma}, where \code{sigma} is the estimated residual standard deviation.
The y-axis therefore corresponds to quarter-root variances.
The y-axis was changed from log2-variance to quarter-root variance in limma version 3.31.21.
The quarter-root scale matches the similar plot produced by the \code{voom} function and gives a better plot when some of the variances are close to zero.

See \code{\link[graphics]{points}} for possible values for \code{pch} and \code{cex}.
}

\value{A plot is created on the current graphics device.}
\author{Gordon Smyth}
\seealso{
\code{\link{eBayes}}

An overview of diagnostic functions available in LIMMA is given in \link{09.Diagnostics}.
}
\keyword{hplot}
