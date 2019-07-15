\name{topRomer}
\alias{topRomer}
\title{Top Gene Set Testing Results from Romer}
\description{
Extract a matrix of the top gene set testing results from the \link{romer} output.
}
\usage{
topRomer(x,n=10,alternative="up")
}
\arguments{
  \item{x}{matrix which is the output from \link{romer}.}
  \item{n}{number of top gene set testing results to be extracted.}
  \item{alternative}{character which can be one of the three possible alternative p values: "up", "down" or "mixed". }
}
\value{
matrix, which is sorted by the "up", "down" or "mixed" p values, with the rows corresponding to estimated p-values for the top number of gene sets and 
the columns corresponding to the number of genes for each gene set and the alternative hypotheses mixed, up, down.
}
\details{
This function takes the results from romer and returns a number of top gene set testing results that are sorted by the p values.
}

\examples{
# See romer for examples
}

\seealso{
\code{\link{romer}}

There is a topic page on \link{10.GeneSetTests}.
}

\author{Gordon Smyth and Yifang Hu}

\keyword{gene set test}
