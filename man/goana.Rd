\name{goana}
\alias{goana}
\alias{goana.default}
\alias{goana.MArrayLM}
\alias{kegga}
\alias{kegga.default}
\alias{kegga.MArrayLM}
\alias{getGeneKEGGLinks}
\alias{getKEGGPathwayNames}
\title{Gene Ontology or KEGG Pathway Analysis}

\description{
Test for over-representation of gene ontology (GO) terms or KEGG pathways in one or more sets of genes, optionally adjusting for abundance or gene length bias.
}

\usage{
\method{goana}{MArrayLM}(de, coef = ncol(de), geneid = rownames(de), FDR = 0.05, trend = FALSE, \dots)
\method{kegga}{MArrayLM}(de, coef = ncol(de), geneid = rownames(de), FDR = 0.05, trend = FALSE, \dots)
\method{goana}{default}(de, universe = NULL, species = "Hs", prior.prob = NULL, covariate=NULL,
      plot=FALSE, \dots)
\method{kegga}{default}(de, universe = NULL, restrict.universe = FALSE, species = "Hs", species.KEGG = NULL,
      convert = FALSE, gene.pathway = NULL, pathway.names = NULL,
      prior.prob = NULL, covariate=NULL, plot=FALSE, \dots)
getGeneKEGGLinks(species.KEGG = "hsa", convert = FALSE)
getKEGGPathwayNames(species.KEGG = NULL, remove.qualifier = FALSE)
}

\arguments{
  \item{de}{a character vector of Entrez Gene IDs, or a list of such vectors, or an \code{MArrayLM} fit object.}

  \item{coef}{column number or column name specifying for which coefficient or contrast differential expression should be assessed.}

  \item{geneid}{Entrez Gene identifiers. Either a vector of length \code{nrow(de)} or the name of the column of \code{de$genes} containing the Entrez Gene IDs.}

  \item{FDR}{false discovery rate cutoff for differentially expressed genes. Numeric value between 0 and 1.}

  \item{species}{character string specifying the species.
  Possible values include \code{"Hs"} (human), \code{"Mm"} (mouse), \code{"Rn"} (rat), \code{"Dm"} (fly) or \code{"Pt"} (chimpanzee), but other values are possible if the corresponding organism package is available.
  See \code{\link{alias2Symbol}} for other possible values.
  Ignored if \code{species.KEGG} or is not \code{NULL} or if \code{gene.pathway} and \code{pathway.names} are not \code{NULL}.}

  \item{species.KEGG}{three-letter KEGG species identifier. See \url{http://www.kegg.jp/kegg/catalog/org_list.html} or \url{http://rest.kegg.jp/list/organism} for possible values.
  Alternatively, if \code{de} contains KEGG ortholog Ids (\code{"k00001"} etc) instead of gene Ids, then set \code{species.KEGG="ko"}.
  This argument is ignored if \code{gene.pathway} and \code{pathway.names} are both not \code{NULL}.}

  \item{convert}{if \code{TRUE} then KEGG gene identifiers will be converted to NCBI Entrez Gene identifiers.  Note that KEGG IDs are the same as Entrez Gene IDs for most species anyway.}

  \item{gene.pathway}{data.frame linking genes to pathways.  First column gives gene IDs, second column gives pathway IDs. By default this is obtained automatically by \code{getGeneKEGGLinks(species.KEGG)}.}

  \item{remove.qualifier}{if \code{TRUE}, the species qualifier will be removed from the pathway names.}

  \item{pathway.names}{data.frame giving full names of pathways. First column gives pathway IDs, second column gives pathway names. By default this is obtained automatically using \code{getKEGGPathwayNames(species.KEGG, remove=TRUE)}.}

  \item{trend}{adjust analysis for gene length or abundance?
  Can be logical, or a numeric vector of covariate values, or the name of the column of \code{de$genes} containing the covariate values.
  If \code{TRUE}, then \code{de$Amean} is used as the covariate.}

  \item{universe}{vector specifying the set of Entrez Gene identifiers to be the background universe.
  If \code{NULL} then all Entrez Gene IDs associated with any gene ontology term will be used as the universe.}

  \item{restrict.universe}{logical, should the \code{universe} be restricted to gene identifiers found in at least one pathway in \code{gene.pathway}?}

  \item{prior.prob}{optional numeric vector of the same length as \code{universe} giving the prior probability that each gene in the universe appears in a gene set.
  Will be computed from \code{covariate} if the latter is provided.
  Ignored if \code{universe} is \code{NULL}.}

  \item{covariate}{optional numeric vector of the same length as \code{universe} giving a covariate against which \code{prior.prob} should be computed.
  Ignored if \code{universe} is \code{NULL}.}

  \item{plot}{logical, should the \code{prior.prob} vs \code{covariate} trend be plotted?}

  \item{\dots}{any other arguments in a call to the \code{MArrayLM} methods are passed to the corresponding default method.}
}

\details{
These functions perform over-representation analyses for Gene Ontology terms or KEGG pathways.
The default methods accept a gene set as a vector of Entrez Gene IDs or multiple gene sets as a list of such vectors.
An over-represention analysis is then done for each set.
The \code{MArrayLM} method extracts the gene sets automatically from a linear model fit object.

The p-values returned by \code{goana} and \code{kegga} are unadjusted for multiple testing.
The authors have chosen not to correct automatically for multiple testing because GO terms and KEGG pathways are often overlapping, so standard methods of p-value adjustment may be very conservative.
Users should be aware though that p-values are unadjusted, meaning that only very small p-values should be used for published results.

\code{goana} uses annotation from the appropriate Bioconductor organism package.
The \code{species} can be any character string XX for which an organism package org.XX.eg.db is installed.
Examples are \code{"Hs"} for human for {"Mm"} for mouse.
See \code{\link{alias2Symbol}} for other possible values for \code{species}.

\code{kegga} reads KEGG pathway annotation from the KEGG website.
For \code{kegga}, the species name can be provided in either Bioconductor or KEGG format.
Examples of KEGG format are \code{"hsa"} for human, \code{"mmu"} for mouse of \code{"dme"} for fly.
\code{kegga} can be used for any species supported by KEGG, of which there are more than 14,000 possibilities.
By default, \code{kegga} obtains the KEGG annotation for the specified species from the \url{http://rest.kegg.jp} website.
Alternatively one can supply the required pathway annotation to \code{kegga} in the form of two data.frames.
If this is done, then an internet connection is not required.

The gene ID system used by \code{kegga} for each species is determined by KEGG.
For human and mouse, the default (and only choice) is Entrez Gene ID.
For Drosophila, the default is FlyBase CG annotation symbol.
The format of the IDs can be seen by typing \code{head(getGeneKEGGLinks(species))}, for example\code{head(getGeneKEGGLinks("hsa"))} or \code{head(getGeneKEGGLinks("dme"))}.
Entrez Gene IDs can always be used.
If Entrez Gene IDs are not the default, then conversion can be done by specifying \code{"convert=TRUE"}.

Another possibility is to use KEGG orthology IDs as the gene IDs, and these can be used for any species.
In that case, set \code{species.KEGG="ko"}.

The ability to supply data.frame annotation to \code{kegga} means that \code{kegga} can in principle be used in conjunction with any user-supplied set of annotation terms.

The default \code{goana} and \code{kegga} methods accept a vector \code{prior.prob} giving the prior probability that each gene in the universe appears in a gene set.
This vector can be used to correct for unwanted trends in the differential expression analysis associated with gene length, gene abundance or any other covariate (Young et al, 2010).
The \code{MArrayLM} object computes the \code{prior.prob} vector automatically when \code{trend} is non-\code{NULL}.

If \code{prior.prob=NULL}, the function computes one-sided hypergeometric tests equivalent to Fisher's exact test.
If prior probabilities are specified, then a test based on the Wallenius' noncentral hypergeometric distribution is used to adjust for the relative probability that each gene will appear in a gene set, following the approach of Young et al (2010).

The \code{MArrayLM} methods performs over-representation analyses for the up and down differentially expressed genes from a linear model analysis.
In this case, the universe is all the genes found in the fit object.

\code{trend=FALSE} is equivalent to \code{prior.prob=NULL}.
If \code{trend=TRUE} or a covariate is supplied, then a trend is fitted to the differential expression results and this is used to set \code{prior.prob}.

The statistical approach provided here is the same as that provided by the goseq package, with one methodological difference and a few restrictions.
Unlike the goseq package, the gene identifiers here must be Entrez Gene IDs and the user is assumed to be able to supply gene lengths if necessary.
The goseq package has additional functionality to convert gene identifiers and to provide gene lengths.
The only methodological difference is that \code{goana} and \code{kegga} computes gene length or abundance bias using \code{\link{tricubeMovingAverage}} instead of monotonic regression.
While \code{\link{tricubeMovingAverage}} does not enforce monotonicity, it has the advantage of numerical stability when \code{de} contains only a small number of genes.
}

\note{
\code{kegga} requires an internet connection unless \code{gene.pathway} and \code{pathway.names} are both supplied.

The default for \code{kegga} with \code{species="Dm"} changed from \code{convert=TRUE} to \code{convert=FALSE} in limma 3.27.8.
Users wanting to use Entrez Gene IDs for Drosophila should set \code{convert=TRUE}, otherwise fly-base CG annotation symbol IDs are assumed (for example "Dme1_CG4637").

The default for \code{restrict.universe=TRUE} in \code{kegga} changed from \code{TRUE} to \code{FALSE} in limma 3.33.4.

Bug fix: results from \code{kegga} with \code{trend=TRUE} or with non-NULL \code{covariate} were incorrect prior to limma 3.32.3.
The results were biased towards significant Down p-values and against significant Up p-values.
}

\value{
The \code{goana} default method produces a data frame with a row for each GO term and the following columns:
  \item{Term}{GO term.}
  \item{Ont}{ontology that the GO term belongs to.  Possible values are \code{"BP"}, \code{"CC"} and \code{"MF"}.}
  \item{N}{number of genes in the GO term.}
  \item{DE}{number of genes in the \code{DE} set.}
  \item{P.DE}{p-value for over-representation of the GO term in the set.}
The last two column names above assume one gene set with the name \code{DE}.
In general, there will be a pair of such columns for each gene set and the name of the set will appear in place of \code{"DE"}.

The \code{goana} method for \code{MArrayLM} objects produces a data frame with a row for each GO term and the following columns:
  \item{Term}{GO term.}
  \item{Ont}{ontology that the GO term belongs to.  Possible values are \code{"BP"}, \code{"CC"} and \code{"MF"}.}
  \item{N}{number of genes in the GO term.}
  \item{Up}{number of up-regulated differentially expressed genes.}
  \item{Down}{number of down-regulated differentially expressed genes.}
  \item{P.Up}{p-value for over-representation of GO term in up-regulated genes.}
  \item{P.Down}{p-value for over-representation of GO term in down-regulated genes.}

The row names of the data frame give the GO term IDs.

The output from \code{kegga} is the same except that row names become KEGG pathway IDs, \code{Term} becomes \code{Pathway} and there is no \code{Ont} column.
}

\references{
  Young, M. D., Wakefield, M. J., Smyth, G. K., Oshlack, A. (2010).
  Gene ontology analysis for RNA-seq: accounting for selection bias.
  \emph{Genome Biology} 11, R14.
  \url{http://genomebiology.com/2010/11/2/R14}
}

\seealso{
\code{\link{topGO}}, \code{\link{topKEGG}}

The goseq package provides an alternative implementation of methods from Young et al (2010).
Unlike the limma functions documented here, goseq will work with a variety of gene identifiers and includes a database of gene length information for various species.

The gostats package also does GO analyses without adjustment for bias but with some other options.

See \link{10.GeneSetTests} for a description of other functions used for gene set testing.
}

\author{Gordon Smyth and Yifang Hu}

\examples{
\dontrun{
## Linear model usage:

fit <- lmFit(y, design)
fit <- eBayes(fit)

# Standard GO analysis

go.fisher <- goana(fit, species="Hs")
topGO(go.fisher, sort = "up")
topGO(go.fisher, sort = "down")

# GO analysis adjusting for gene abundance

go.abund <- goana(fit, geneid = "GeneID", trend = TRUE)
topGO(go.abund, sort = "up")
topGO(go.abund, sort = "down")

# GO analysis adjusting for gene length bias
# (assuming that y$genes$Length contains gene lengths)

go.len <- goana(fit, geneid = "GeneID", trend = "Length")
topGO(go.len, sort = "up")
topGO(go.len, sort = "down")

## Default usage with a list of gene sets:

go.de <- goana(list(DE1 = EG.DE1, DE2 = EG.DE2, DE3 = EG.DE3))
topGO(go.de, sort = "DE1")
topGO(go.de, sort = "DE2")
topGO(go.de, ontology = "BP", sort = "DE3")
topGO(go.de, ontology = "CC", sort = "DE3")
topGO(go.de, ontology = "MF", sort = "DE3")

## Standard KEGG analysis

k <- kegga(fit, species="Hs")
k <- kegga(fit, species.KEGG="hsa") # equivalent to previous
topKEGG(k, sort = "up")
topKEGG(k, sort = "down")

}
}

\keyword{gene set test}
