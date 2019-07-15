\name{cumOverlap}
\alias{cumOverlap}
\title{Cumulative Overlap Analysis of Ordered Lists}
\description{
Test whether the leading members of ordered lists significantly overlap.
}
\usage{
cumOverlap(ol1, ol2)
}

\arguments{
  \item{ol1}{vector containing first ordered list. Duplicate values not allowed.}
  \item{ol2}{vector containing second ordered list. Should contain the same values as found in \code{ol1} but in a possibly different order. Duplicate values not allowed.}
}

\value{
List containing the following components:
\item{n.total}{integer, total number of values in common between \code{ol1} and \code{ol2}.}
\item{n.min}{integer, top table length leading to smallest adjusted p-value.}
\item{p.min}{smallest adjusted p-value.}
\item{n.overlap}{integer, number of overlapping IDs in first \code{n.min}.}
\item{id.overlap}{vector giving the overlapping IDs in first \code{n.min}.}
\item{p.value}{numeric, vector of p-values for each possible top table length.}
\item{adj.p.value}{numeric, vector of Bonferroni adjusted p-values for each possible top table length.}
}

\details{
The function compares the top \code{n} members of each list, for every possible \code{n}, and conducts an hypergeometric test for overlap.
The function returns the value of \code{n} giving the smallest p-value.

The p-values are adjusted for multiple testing in a similar way to Bonferroni's method, but starting from the top of th e ranked list instead of from the smallest p-values.
This approach is designed to be sensitive to contexts where the number of Ids involved in the significant overlap are a small proportion of the total.

The vectors \code{ol1} and \code{ol2} do not need to be of the same length, but only values in common between the two vectors will be used in the calculation.

This method was described in Chapter 4 of Wu (2011).
}

\author{Gordon Smyth and Di Wu}

\references{
Wu, D (2011).
Finding hidden relationships between gene expression profiles with application to breast cancer biology.
PhD thesis, University of Melbourne.
\url{http://hdl.handle.net/11343/36278}
}

\examples{
GeneIds <- paste0("Gene",1:50)
ol1 <- GeneIds
ol2 <- c(sample(GeneIds[1:5]), sample(GeneIds[6:50]))
coa <- cumOverlap(ol1, ol2)
coa$p.min
coa$id.overlap
}

\concept{gene set enrichment analysis}
