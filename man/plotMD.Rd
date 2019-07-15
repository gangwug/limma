\title{Mean-Difference Plot of Expression Data}
\name{plotMD}
\alias{plotMD}
\alias{plotMD.RGList}
\alias{plotMD.MAList}
\alias{plotMD.EListRaw}
\alias{plotMD.EList}
\alias{plotMD.MArrayLM}
\alias{plotMD.default}

\description{
Creates a mean-difference plot (aka MA plot) with color coding for highlighted points.
}

\usage{
\method{plotMD}{default}(object, column = 1, xlab = "Average log-expression",
       ylab = "Expression log-ratio (this sample vs others)",
       main = colnames(object)[column], status=NULL, \dots)
\method{plotMD}{EList}(object, column = 1, array = NULL, xlab = "Average log-expression",
       ylab = "Expression log-ratio (this sample vs others)",
       main = colnames(object)[column], status=object$genes$Status,
       zero.weights = FALSE, \dots)
\method{plotMD}{RGList}(object, column = 1, array = NULL, xlab = "A", ylab = "M",
       main = colnames(object)[column], status=object$genes$Status,
       zero.weights = FALSE, \dots)
\method{plotMD}{MAList}(object, column = 1, array = NULL, xlab = "A", ylab = "M",
       main = colnames(object)[column], status=object$genes$Status,
       zero.weights = FALSE, \dots)
\method{plotMD}{MArrayLM}(object, column = ncol(object), coef = NULL, xlab = "Average log-expression",
       ylab = "log-fold-change", main = colnames(object)[column],
       status=object$genes$Status, zero.weights = FALSE, \dots)
}

\arguments{
  \item{object}{an \code{RGList}, \code{MAList}, \code{EList}, \code{ExpressionSet} or \code{MArrayLM} object.
  Alternatively a numeric \code{matrix}.}
  \item{column}{integer, column of \code{object} to be plotted.}
  \item{array}{alternative to \code{column} for microarray data objects. If specified, then \code{column} is ignored.}
  \item{coef}{alternative to \code{column} for fitted model objects. If specified, then \code{column} is ignored.}
  \item{xlab}{character string, label for x-axis}
  \item{ylab}{character string, label for y-axis}
  \item{main}{character string, title for plot}
  \item{status}{vector giving the control status of each spot on the array, of same length as the number of rows of \code{object}.
  If \code{NULL}, then all points are plotted in the default color, symbol and size.}
  \item{zero.weights}{logical, should spots with zero or negative weights be plotted?}
  \item{\dots}{other arguments are passed to \code{\link{plotWithHighlights}}.}
}

\details{
A mean-difference plot (MD-plot) is a plot of log-intensity ratios (differences) versus log-intensity averages (means).
For two color data objects, a within-array MD-plot is produced with the M and A values computed from the two channels for the specified array.
This is the same as a mean-difference plot (\code{\link{mdplot}}) with the red and green log2-intensities of the array providing the two columns.

For single channel data objects, a between-array MD-plot is produced.
An articifial array is produced by averaging all the arrays other than the array specified.
A mean-difference plot is then producing from the specified array and the artificial array.
Note that this procedure reduces to an ordinary mean-difference plot when there are just two arrays total.

If \code{object} is an \code{MArrayLM} object, then the plot is an fitted model MD-plot in which the estimated coefficient is on the y-axis and the average A-value is on the x-axis.

The \code{status} vector can correspond to any grouping of the probes that is of interest.
If \code{object} is a fitted model object, then \code{status} vector is often used to indicate statistically significance, so that differentially expressed points are highlighted.
If \code{object} is a microarray data object, then \code{status} might distinguish control probes from regular probes so that different types of controls are highlighted.

The \code{status} can be included as the component \code{object$genes$Status} instead of being passed as an argument to \code{plotMD}.

See \code{\link{plotWithHighlights}} for how to set colors and graphics parameters for the highlighted and non-highlighted points.
}

\note{
This function is an alternative to \code{plotMA}, which was one of the original functions of the limma package in 2002.
The history of mean-difference plots and MA-plots is reviewed in Ritchie et al (2015).
}

\value{A plot is created on the current graphics device.}

\references{
Ritchie, ME, Phipson, B, Wu, D, Hu, Y, Law, CW, Shi, W, and Smyth, GK (2015).
limma powers differential expression analyses for RNA-sequencing and microarray studies.
\emph{Nucleic Acids Research} Volume 43, e47.
\url{http://nar.oxfordjournals.org/content/43/7/e47}
}

\author{Gordon Smyth}

\examples{
A <- runif(1000,4,16)
y <- A + matrix(rnorm(1000*3,sd=0.2),1000,3)
status <- rep(c(0,-1,1),c(950,40,10))
y[,1] <- y[,1] + status
plotMD(y, column=1, status=status, values=c(-1,1), hl.col=c("blue","red"))

MA <- new("MAList")
MA$A <- runif(300,4,16)
MA$M <- rt(300,df=3)

# Spike-in values
MA$M[1:3] <- 0
MA$M[4:6] <- 3
MA$M[7:9] <- -3

status <- rep("Gene",300)
status[1:3] <- "M=0"
status[4:6] <- "M=3"
status[7:9] <- "M=-3"
values <- c("M=0","M=3","M=-3")
hl.col <- c("blue","red","green")

plotMD(MA,main="MA-Plot with 12 spiked-in points",
       status=status, values=values, hl.col=hl.col)

#  Same as above but setting graphical parameters as attributes
attr(status,"values") <- values
attr(status,"col") <- hl.col
plotMD(MA, main="Mean-Difference Plot with 12 spiked-in points", status=status)

#  Same as above but passing status as part of object
MA$genes$Status <- status
plotMD(MA, main="Mean-Difference Plot with 12 spiked-in points")

#  Change settings for background points
MA$genes$Status <- status
plotMD(MA, bg.pch=1, bg.cex=0.5)
}

\seealso{
The driver function for \code{plotMD} is \code{\link{plotWithHighlights}}.
See also \code{\link{mdplot}} for a very basic mean-difference plot function.

An overview of plot functions available in LIMMA is given in \link{09.Diagnostics}.
}

\keyword{hplot}
