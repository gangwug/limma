\name{topGO}
\alias{topGO}
\alias{topKEGG}
\title{Table of Top GO Terms or Top KEGG Pathways}

\description{
Extract top GO terms from goana output or top KEGG pathways from kegga output.
}

\usage{
topGO(results, ontology = c("BP", "CC", "MF"), sort = NULL, number = 20L, 
      truncate.term = NULL)
topKEGG(results, sort = NULL, number = 20L, truncate.path = NULL)
}

\arguments{
  \item{results}{data frame produced by \code{\link{goana}} or \code{\link{kegga}}.} 
  \item{ontology}{character vector of ontologies to be included in output.  Elements should be one or more of \code{"BP"}, \code{"CC"} or \code{"MF"}.}
  \item{sort}{character vector of names of gene lists for which results are required.  Should be one or more of the column names of \code{results}.  Defaults to all gene lists.}
  \item{number}{maximum number of top GO terms or top KEGG pathways to list. For all terms or all pathways, set \code{number=Inf}.}
  \item{truncate.term}{truncate the name of the GO term at this number of characters.}
  \item{truncate.path}{truncate the name of the KEGG pathway at this number of characters.}
}

\details{
\code{topGO} organizes the output from \code{\link{goana}} into top-tables of the most significant GO terms.
\code{topKEGG} similarly extracts the most significant KEGG pathways from \code{\link{kegga}} output. 
In either case, rows are sorted by the minimum p-value of any of the result columns specified by \code{sort}.
}

\value{
Same as \code{results} but with rows subsetted by Ontology and sorted by p-value.
}

\seealso{
\code{\link{goana}}, \code{\link{kegga}}

See \link{10.GeneSetTests} for a description of other functions used for gene set testing.
}

\author{Gordon Smyth and Yifang Hu}

\examples{
# See goana examples
}

\keyword{gene set test}
