\name{barcodeplot}
\alias{barcodeplot}
\title{Barcode Enrichment Plot}
\description{
Display the enrichment of one or two gene sets in a ranked gene list.
}
\usage{
barcodeplot(statistics, index = NULL, index2 = NULL, gene.weights = NULL,
            weights.label = "Weight", labels = c("Down","Up"),
            quantiles = c(-1,1)*sqrt(2), col.bars = NULL, alpha = 0.4,
            worm = TRUE, span.worm = 0.45, xlab = "Statistic", \dots)
}
\arguments{
  \item{statistics}{numeric vector giving the values of statistics to rank genes by.}
  \item{index}{index vector for the gene set.
  This can be a vector of indices, or a logical vector of the same length as \code{statistics} or, in general, any vector such that \code{statistic[index]} gives a subset of the statistic values.
  Can be omitted if \code{gene.weights} has same length as \code{statistics}, in which case positive values of \code{gene.weights} indicate to members of the positive set and negative weights correspond to members of the negative set.}
  \item{index2}{optional index vector for a second (negative) gene set.
  If specified, then \code{index} and \code{index2} specify positive and negative genes respectively.
  Usually used to distinguish down-regulated genes from up-regulated genes.}
  \item{gene.weights}{numeric vector giving directional weights for the genes in the (first) set.
  Positive and negative weights correspond to positive and negative genes.
  Ignored if \code{index2} is non-null.}
  \item{weights.label}{label describing the entries in \code{gene.weights}.}
  \item{labels}{character vector of labels for low and high statistics.  First label is associated with low statistics or negative statistics and is displayed at the left end of the plot.  Second label is associated with high or positive statistics and is displayed at the right end of the plot.}
  \item{quantiles}{numeric vector of length 2, giving cutoff values for \code{statistics} considered small or large respectively.  Used to color the rectangle of the barcodeplot.}
  \item{col.bars}{character vector of colors for the vertical bars of the barcodeplot showing the ranks of the gene set members.
  Defaults to \code{"black"} for one set or \code{c("red","blue")} for two sets.}
  \item{alpha}{transparency for vertical bars.  When \code{gene.weights} are not \code{NULL}, values \code{0<alpha<1} give semitransparent colors for the vertical bars inside the rectangle.  This helps distinguish position bars from the weighted bars and also helps to show the density of the bars when there are many bars.  Ignored if \code{gene.weights=NULL}.}
  \item{worm}{logical, should enrichment worms be plotted?}
  \item{span.worm}{loess span for enrichment worms.  Larger spans give smoother worms.}
  \item{xlab}{x-axis label for \code{statistics}.}
  \item{\ldots}{other arguments are passed to \code{plot}.}
}

\value{
No value is returned but a plot is produced as a side effect.
}

\details{
This function plots the positions of one or two gene sets in a ranked list of statistics.
If there are two sets, then one is considered to be the positive set and the other the down set.
For example, the first set and second sets often correspond to genes that are expected to be up- or down-regulated respectively.
The function can optionally display varying weights for different genes, for example log-fold-changes from a previous experiment.

The statistics are ranked left to right from smallest to largest.
The ranked statistics are represented by a shaded bar or bed, and the positions of the specified subsets are marked by vertical bars, forming a pattern like a barcode.
An enrichment worm optionally shows the relative enrichment of the vertical bars in each part of the plot.
The worm is computed by the \code{tricubeMovingAverage} function.

Barcode plots are often used in conjunction with gene set tests, and show the enrichment of gene sets amongst high or low ranked genes.
They were inspired by the set location plot of Subramanian et al (2005), with a number of enhancements, especially the ability to plot positive and negative sets simultaneously.
Barcode plots first appeared in the literature in Lim et al (2009).
More recent examples can be seen in Liu et al (2014), Sheikh et al (2015), Witkowski et al (2015) and Ng et al (2015).

The function can be used with any of four different calling sequences:
  \itemize{
    \item \code{index} is specified, but not \code{index2} or \code{gene.weights}.  Single direction plot.
    \item \code{index} and \code{index2} are specified.  Two directional plot.
    \item \code{index} and \code{gene.weights} are specified.  \code{gene.weights} must have same length as \code{statistics[index]}.  Plot will be two-directional if \code{gene.weights} contains positive and negative values.
    \item \code{gene.weights} is specified by not \code{index} or \code{index2}.  \code{gene.weights} must have same length as \code{statistics}.  Plot will be two-directional if \code{gene.weights} contains positive and negative values.      
  }
}

\seealso{
\code{\link{tricubeMovingAverage}}, \code{\link{roast}}, \code{\link{camera}}, \code{\link{romer}}, \code{\link{geneSetTest}}

There is a topic page on \link{10.GeneSetTests}.
}

\author{Yifang Hu, Gordon Smyth and Di Wu}

\references{
Ng, AP, Hu, Y, Metcalf, D, Hyland, CD, Ierino, H, Phipson, B, Wu, D, Baldwin, TM, Kauppi, M, Kiu, H, Di, Rago, L, Hilton, DJ, Smyth, GK, Alexander, WS (2015).
Early lineage priming by trisomy of Erg leads to myeloproliferation in a down syndrome model.
\emph{PLOS Genetics} 11, e1005211. 
\url{http://www.ncbi.nlm.nih.gov/pubmed/25973911}

Lim E, Vaillant F, Wu D, Forrest NC, Pal B, Hart AH, Asselin-Labat ML, Gyorki DE, Ward T, Partanen A, Feleppa F, Huschtscha LI, Thorne HJ; kConFab; Fox SB, Yan M, French JD, Brown MA, Smyth GK, Visvader JE, and Lindeman GJ (2009).
Aberrant luminal progenitors as the candidate target population for basal tumor development in BRCA1 mutation carriers.
\emph{Nature Medicine} 15, 907-913.

Liu, GJ, Cimmino, L, Jude, JG, Hu, Y, Witkowski, MT, McKenzie, MD, Kartal-Kaess, M, Best, SA, Tuohey, L, Liao, Y, Shi, W, Mullighan, CG, Farrar, MA, Nutt, SL, Smyth, GK, Zuber, J, and Dickins, RA (2014).
Pax5 loss imposes a reversible differentiation block in B progenitor acute lymphoblastic leukemia.
\emph{Genes & Development} 28, 1337-1350. 
\url{http://www.ncbi.nlm.nih.gov/pubmed/24939936}

Sheikh, B, Lee, S, El-saafin, F, Vanyai, H, Hu, Y, Pang, SHM, Grabow, S, Strasser, A, Nutt, SL, Alexander, WS, Smyth, GK, Voss, AK, and Thomas, T (2015). 
MOZ regulates B cell progenitors in mice, consequently, Moz haploinsufficiency dramatically retards MYC-induced lymphoma development. 
\emph{Blood} 125, 1910-1921.
\url{http://www.ncbi.nlm.nih.gov/pubmed/25605372}

Subramanian A, Tamayo P, Mootha VK, Mukherjee S, Ebert BL, Gillette MA, Paulovich A, Pomeroy SL, Golub TR, Lander ES, and Mesirov JP (2005).
Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles.
\emph{Proc Natl Acad Sci USA} 102, 15545-15550.

Witkowski, MT, Cimmino, L, Hu, Y, Trimarchi, T, Tagoh, H, McKenzie, MD, Best, SA, Tuohey, L, Willson, TA, Nutt, SL, Meinrad Busslinger, M, Aifantis, I, Smyth, GK, and Dickins, RA (2015). 
Activated Notch counteracts Ikaros tumor suppression in mouse and human T cell acute lymphoblastic leukemia.
\emph{Leukemia} 29, 1301-1311.
\url{http://www.ncbi.nlm.nih.gov/pubmed/25655195}
}

\examples{
stat <- rnorm(100)
sel <- 1:10
sel2 <- 11:20
stat[sel] <- stat[sel]+1
stat[sel2] <- stat[sel2]-1

# One directional
barcodeplot(stat, index = sel)

# Two directional
barcodeplot(stat, index = sel, index2 = sel2)

# Second set can be indicated by negative weights
barcodeplot(stat, index = c(sel,sel2), gene.weights = c(rep(1,10), rep(-1,10)))

# Two directional with unequal weights
w <- rep(0,100)
w[sel] <- runif(10)
w[sel2] <- -runif(10)
barcodeplot(stat, gene.weights = w, weights.label = "logFC")

# One directional with unequal weights
w <- rep(0,100)
w[sel2] <- -runif(10)
barcodeplot(stat, gene.weights = w, weights.label = "logFC", col.bars = "dodgerblue")
}
