\title{Heatmap of gene expression values}
\name{coolmap}
\alias{coolmap}

\description{
Create a heatmap of a matrix of log-expression values.
}

\usage{
coolmap(x, cluster.by="de pattern", col=NULL,
        linkage.row="complete", linkage.col="complete", show.dendrogram="both", ...)
}

\arguments{
  \item{x}{any data object that can be coerced to a matrix of log-expression values, for example an \code{ExpressionSet} or \code{EList}. Rows represent genes and columns represent RNA samples.}
  \item{cluster.by}{choices are \code{"de pattern"} or \code{"expression level"}.
  In the former case, the intention is to cluster by relative changes in expression, so genes are clustered by Pearson correlation and log-expression values are mean-corrected by rows for the plot.
  In the latter case, the intention is to cluster by absolute expression, so genes are clustered by Euclidean and log-expression values are not mean-corrected.}
  \item{col}{character vector specifying the color panel.
  Can be either the name of the panel or a vector of R colors that can be passed directly to the \code{heatmap.2} function.
  Possible panel names are \code{"redblue"}, \code{"redgreen"}, \code{"yellowblue"} or \code{"whitered"}.
  Defaults to \code{"redblue"} if \code{cluster.by="de pattern"} or \code{"yellowblue"} if \code{cluster.by="expression level"}.}
  \item{linkage.row}{linkage criterion used to cluster the rows.
  Choices are \code{"none"}, \code{"ward"}, \code{"ward.D"}, \code{"ward.D2"}, \code{"single"}, \code{"complete"}, \code{"average"}, \code{"mcquitty"}, \code{"median"} or \code{"centroid"}, with \code{"ward"} treated as \code{"ward.D2"}.}
  \item{linkage.col}{linkage criterion used to cluster the columns.
  Choices are the same as for \code{linkage.row}.}
  \item{show.dendrogram}{choices are \code{"row"}, \code{"column"}, \code{"both"} or \code{"none"}.}
  \item{\dots}{any other arguments are passed to \code{heatmap.2}.}
}

\details{
This function calls the \code{heatmap.2} function in the ggplots package with sensible argument settings for genomic log-expression data.
The default settings for \code{heatmap.2} are often not ideal for expression data, and overriding the defaults requires explicit calls to \code{hclust} and \code{as.dendrogram} as well as prior standardization of the data values.
The \code{coolmap} function implements our preferred defaults for the two most common types of heatmaps.
When clustering by relative expression (\code{cluster.by="de pattern"}), it implements a row standardization that takes account of \code{NA} values and standard deviations that might be zero.
}

\value{
A plot is created on the current graphics device.
A list is also invisibly returned, see \code{\link[gplots]{heatmap.2}} for details.
}

\author{Gordon Smyth}

\seealso{
\code{\link[gplots]{heatmap.2}}, \code{\link{hclust}}, \code{\link{dist}}.

An overview of diagnostic functions available in LIMMA is given in \link{09.Diagnostics}.
}

\examples{
# Simulate gene expression data for 50 genes and 6 microarrays.
# Samples are in two groups
# First 50 probes are differentially expressed in second group
ngenes <- 50
sd <- 0.3*sqrt(4/rchisq(ngenes,df=4))
x <- matrix(rnorm(ngenes*6,sd=sd),ngenes,6)
rownames(x) <- paste("Gene",1:ngenes)
x <- x + seq(from=0, to=16, length=ngenes)
x[,4:6] <- x[,4:6] + 2
coolmap(x)
}
