\name{ids2indices}
\alias{ids2indices}
\title{Convert Gene Identifiers to Indices for Gene Sets}
\description{
Make a list of gene identifiers into a list of indices for gene sets.
}
\usage{
ids2indices(gene.sets, identifiers, remove.empty=TRUE)
}
\arguments{
  \item{gene.sets}{list of character vectors, each vector containing the gene identifiers for a set of genes.}
  \item{identifiers}{character vector of gene identifiers.}
  \item{remove.empty}{logical, should sets of size zero be removed from the output?}
}
\value{
list of integer vectors, each vector containing the indices of a gene set in the vector \code{identifiers}.
}
\details{
This function used to create input for \code{romer}, \code{mroast} and \code{camera} function.
Typically, \code{identifiers} is the vector of Entrez Gene IDs, and \code{gene.sets} is obtained constructed from a database of gene sets, 
for example a representation of the Molecular Signatures Database (MSigDB) downloaded from \url{http://bioinf.wehi.edu.au/software/MSigDB}.
}

\seealso{
\code{\link{romer}}, \code{\link{mroast}}, \code{\link{camera}}

There is a topic page on \link{10.GeneSetTests}.
}

\examples{

\dontrun{

download.file("http://bioinf.wehi.edu.au/software/MSigDB/human_c2_v5p2.rdata", 
	      "human_c2_v5p2.rdata", mode = "wb")

load("human_c2_v5p2.rdata")
c2.indices <- ids2indices(Hs.c2, y$genes$GeneID)
camera(y, c2.indices, design)

}

}
\author{Gordon Smyth and Yifang Hu}

\keyword{gene set test}
