\name{wsva}
\alias{wsva}
\title{Weighted Surrogate Variable Analysis}
\description{
Calculate surrogate variables from the singular vectors of the linear model residual space.
}

\usage{
wsva(y, design, n.sv = 1L, weight.by.sd = FALSE, plot = FALSE, ...)
}

\arguments{
  \item{y}{numeric matrix giving log-expression or log-ratio values for a series of microarrays, or any object that can coerced to a matrix including \code{ExpressionSet}, \code{MAList}, \code{EList} or \code{PLMSet} objects.
        Rows correspond to genes and columns to samples.}
  \item{design}{design matrix}
  \item{n.sv}{number of surrogate variables required.}
  \item{weight.by.sd}{logical, should the surrogate variables be especially tuned to the more variable genes?}
  \item{plot}{logical. If \code{TRUE}, plots the proportion of variance explained by each surrogate variable.}
  \item{\dots}{other arguments can be included that would be suitable for \code{lmFit}.}
}

\value{
Numeric matrix with \code{ncol(y)} rows and \code{n.sv} columns containing the surrogate variables.
}

\details{
The function constructs surrogate variables that explain a high proportion of the residual variability for many of the genes.
The surrogate variables can be included in the design matrix to remove unwanted variation.
The surrogate variables are constructed from the singular vectors of a representation of the linear model residual space.

If \code{weight.by.sd=FALSE}, then the method is a simplification of the approach by Leek and Storey (2007).
}

\author{Gordon Smyth and Yifang Hu}

\references{
Leek, JT, Storey, JD (2007).
Capturing heterogeneity in gene expression studies by surrogate variable analysis.
\emph{PLoS Genetics} 3, 1724-1735.
}
