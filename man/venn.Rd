\name{venn}
\alias{vennCounts}
\alias{vennDiagram}
\title{Venn Diagrams}
\description{
Compute classification counts and draw a Venn diagram.
}
\usage{
vennCounts(x, include="both")
vennDiagram(object, include="both", names=NULL, mar=rep(1,4), cex=c(1.5,1,0.7), lwd=1,
            circle.col=NULL, counts.col=NULL, show.include=NULL, \dots)
}
\arguments{
  \item{x}{a \code{TestResults} matrix.
  This is numeric matrix of 0's, 1's and -1's indicating significance of a test or membership of a set.
  Each row corresponds to a gene and each column to a contrast or set.
  Usually created by \code{\link{decideTests}}.}
  \item{object}{either a \code{TestResults} matrix or a \code{VennCounts} object produced by \code{vennCounts}.}
  \item{include}{character vector specifying whether all differentially expressed genes should be counted, or whether the counts should be restricted to genes changing in a certain direction.  Choices are \code{"both"} for all differentially expressed genes, \code{"up"} for up-regulated genes only or \code{"down"} for down-regulated genes only. If \code{include=c("up","down")} then both the up and down counts will be shown. This argument is ignored if \code{object} if \code{object} is already a \code{vennCounts} object.}
  \item{names}{character vector giving names for the sets or contrasts}
  \item{mar}{numeric vector of length 4 specifying the width of the margins around the plot. This argument is passed to \code{par}.}
  \item{cex}{numerical vector of length 3 giving scaling factors for large, medium and small text on the plot.}
  \item{lwd}{numerical value giving the amount by which the circles should be scaled on the plot. See \code{par}.}
  \item{circle.col}{vector of colors for the circles. See \code{par} for possible values.}
  \item{counts.col}{vector of colors for the counts. Of same length as \code{include}. See \code{par} for possible values.}
  \item{show.include}{logical value whether the value of \code{include} should be printed on the plot. Defaults to \code{FALSE} if \code{include} is a single value and \code{TRUE} otherwise}
  \item{\dots}{any other arguments are passed to \code{plot}}
}

\details{
Each column of \code{x} corresponds to a contrast or set, and the entries of \code{x} indicate membership of each row in each set or alternatively the significance of each row for each contrast.
In the latter case, the entries can be negative as well as positive to indicate the direction of change.

\code{vennCounts} can collate intersection counts for any number of sets.
\code{vennDiagram} can plot up to five sets.
}

\value{
\code{vennCounts} produces an object of class \code{"VennCounts"}.
This contains only one slot, which is numerical matrix with \code{2^ncol{x}} rows and \code{ncol(x)+1} columns.
Each row corresponds to a particular combination of set memberships.
The first \code{ncol{x}} columns of output contain 1 or 0 indicating membership or not in each set.
The last column called \code{"Counts"} gives the number of rows of \code{x} corresponding to that combination of memberships.

\code{vennDiagram} produces no output but causes a plot to be produced on the current graphical device.
}

\seealso{
An overview of linear model functions in limma is given by \link{06.LinearModels}.
}

\author{Gordon Smyth, James Wettenhall, Francois Pepin, Steffen Moeller and Yifang Hu} 

\examples{
Y <- matrix(rnorm(100*6),100,6)
Y[1:10,3:4] <- Y[1:10,3:4]+3
Y[1:20,5:6] <- Y[1:20,5:6]+3
design <- cbind(1,c(0,0,1,1,0,0),c(0,0,0,0,1,1))
fit <- eBayes(lmFit(Y,design))
results <- decideTests(fit)
a <- vennCounts(results)
print(a)
mfrow.old <- par()$mfrow
par(mfrow=c(1,2))
vennDiagram(a)
vennDiagram(results, 
    include=c("up", "down"),
    counts.col=c("red", "blue"),
    circle.col = c("red", "blue", "green3"))
par(mfrow=mfrow.old)
}

\keyword{hplot}
