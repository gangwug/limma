\name{plotWithHighlights}
\alias{plotWithHighlights}
\title{Scatterplot With Highlighting of Special Points}
\description{
Creates scatterplot, with optional size and color coding for points of special interest.
This is the engine for \code{plotMD} and \code{plotMA}.
}
\usage{
plotWithHighlights(x, y, status = NULL, values = NULL,
                   hl.pch = 16, hl.col = NULL, hl.cex = 1, legend = "topright",
                   bg.pch = 16, bg.col = "black", bg.cex = 0.3,
                   pch = NULL, col = NULL, cex = NULL, \dots)
}
\arguments{
  \item{x}{numeric vector.}

  \item{y}{numeric vector.}

  \item{status}{character vector giving the control status of each point, of same length as \code{x} and \code{y}.
  If \code{NULL}, then all points are plotted in the background color, symbol and size.}

  \item{values}{character vector giving values of \code{status} to be highlighted on the plot.
  Defaults to unique values of \code{status} in decreasing order of frequency, with the most frequent value set as the background value.
  Ignored if there is no \code{status} vector.}

  \item{hl.pch}{vector of plotting characters for highlighted points, either of unit length or of same length as \code{values}.
  Ignored is there is no \code{status} vector.}

  \item{hl.col}{vector of colors for highlighted points, either of unit length or of same length as \code{values}.
  Defaults to \code{1+1:length(values)}.
  Ignored if there is no \code{status} vector.}

  \item{hl.cex}{numeric vector of plot symbol expansions for highlighted points, either of unit length or of same length as \code{values}.
  Ignored if there is no \code{status} vector.}

  \item{legend}{character string giving position to place legend.
  See \code{\link{legend}} for possible values.
  Can also be logical, with \code{FALSE} meaning no legend.
  Ignored if there is no \code{status} vector.}

  \item{bg.pch}{plotting character for background (non-highlighted) points.}

  \item{bg.col}{color for background (non-highlighted) points.}

  \item{bg.cex}{plot symbol expansion for background (non-highlighted) points.}

  \item{pch}{synonym for \code{hl.pch} allowed for backward compatibility.}

  \item{col}{synonym for \code{hl.col} allowed for backward compatibility.}

  \item{cex}{synonym for \code{hl.cex} allowed for backward compatibility.}

  \item{\dots}{other arguments are passed to \code{plot}.}
}

\details{
This function produces a scatterplot in which the highlighted points are, by default, larger and colored compared to background points.

The \code{status} vector establishes the status of each point and \code{values} indicates which values of \code{status} should be highlighted.
If \code{values=NULL}, then the most common value of \code{status} is assumed to correspond to background points and all other values are highlighted.

The arguments \code{hl.pch}, \code{hl.col} and \code{hl.cex} give graphics settings for highlighted points.
By default, highlighted points are larger than background points and a different color is used for each distinct highlighted value.

The arguments \code{bg.pch}, \code{bg.col} and \code{bg.cex} give the graphics settings for non-highlighted (background) points.
The same settings are used for all background points.

The arguments \code{values}, \code{pch}, \code{col} and \code{cex} can be included as attributes to \code{status} instead of being passed as arguments to \code{plotWithHighlights}.
This is for compatibility with \code{\link{controlStatus}}.

See \code{\link[graphics]{points}} for possible values for the graphics parameters.
}

\value{A plot is created on the current graphics device.}

\references{
Ritchie, ME, Phipson, B, Wu, D, Hu, Y, Law, CW, Shi, W, and Smyth, GK (2015).
limma powers differential expression analyses for RNA-sequencing and microarray studies.
\emph{Nucleic Acids Research} 43, e47.
\url{http://nar.oxfordjournals.org/content/43/7/e47}
}

\author{Gordon Smyth}

\examples{
x <- runif(1000, min=4, max=16)
status <- rep(c(0,-1,1), c(950,40,10))
y <- status + rnorm(1000, sd=0.2)
plotWithHighlights(x, y, status=status)
}

\seealso{
\code{\link{plotMD}}, \code{\link{plotMA}}, \code{\link{mdplot}}

An overview of diagnostic plots available in LIMMA is given in \link{09.Diagnostics}.
}
