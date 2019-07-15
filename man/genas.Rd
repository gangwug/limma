\name{genas}
\alias{genas}
\alias{genas}
\title{Genuine Association of Gene Expression Profiles}
\description{
Calculates biological correlation between two gene expression profiles.
}
\usage{
genas(fit, coef=c(1,2), subset="all", plot=FALSE, alpha=0.4)
}

\arguments{
 \item{fit}{an \code{MArrayLM} fitted model object produced by \code{lmFit} or \code{contrasts.fit} and followed by \code{eBayes}.}
 \item{coef}{numeric vector of length 2 indicating which columns in the fit object are to be correlated.}
 \item{subset}{character string indicating which subset of genes to include in the correlation analysis.
       Choices are \code{"all"}, \code{"Fpval"}, \code{"p.union"}, \code{"p.int"}, \code{"logFC"} or \code{"predFC"}.}
 \item{plot}{logical, should a scatterplot be produced summarizing the correlation analysis?}
 \item{alpha}{numeric value between 0 and 1 determining the transparency of the technical and biological ellipses if a plot is produced.
 \code{alpha=0} indicates fully transparent and \code{alpha=1} indicates fully opague.}
}

\details{
The function estimates the biological correlation between two different contrasts in a linear model.
By biological correlation, we mean the correlation that would exist between the log2-fold changes (logFC) for the two contrasts, if measurement error could be eliminated and the true log-fold-changes were known.
This function is motivated by the fact that different contrasts for a linear model are often strongly correlated in a technical sense.
For example, the estimated logFC for multiple treatment conditions compared back to the same control group will be positively correlated even in the absence of any biological effect.
This function aims to separate the biological from the technical components of the correlation.
The method is explained briefly in Majewski et al (2010) and in full detail in Phipson (2013).

The \code{subset} argument specifies whether and how the fit object should be subsetted.
Ideally, only genes that are truly differentially expressed for one or both of the contrasts should be used estimate the biological correlation.
The default is \code{"all"}, which uses all genes in the fit object to estimate the biological correlation.
The option \code{"Fpval"} chooses genes based on how many F-test p-values are estimated to be truly significant using the function \code{propTrueNull}.
This should capture genes that display any evidence of differential expression in either of the two contrasts.
The options \code{"p.union"} and \code{"p.int"} are based on the moderated t p-values from both contrasts.
From the \code{propTrueNull} function an estimate of the number of p-values truly significant in either of the two contrasts can be obtained.
"p.union" takes the union of these genes and \code{"p.int"} takes the intersection of these genes.
The other options, \code{"logFC"} and \code{"predFC"} subsets on genes that attain a logFC or predFC at least as large as the 90th percentile of the log fold changes or predictive log fold changes on the absolute scale. 

The \code{plot} option is a logical argument that specifies whether or not to plot a scatter plot of log-fold-changes for the two contrasts.
The biological and technical correlations are overlaid on the scatterplot using semi-transparent ellipses.
\code{library(ellipse)} is required to enable the plotting of ellipses.
}

\note{
As present, \code{genas} assumes that technical correlations between coefficients are the same for all genes, and hence it only works with fit objects that were created without observation weights or missing values.
It does not work with \code{voom} pipelines, because these involve observation weights.
}

\value{
	\code{genas} produces a list with the following components:
	  \item{technical.correlation}{estimate of the technical correlation}
	  \item{biological.correlation}{estimate of the biological correlation}
	  \item{covariance.matrix}{estimate of the covariance matrix from which the biological correlation is obtained}
	  \item{deviance}{the likelihood ratio test statistic used to test whether the biological correlation is equal to 0}
	  \item{p.value}{the p.value associated with \code{deviance}}
	  \item{n}{the number of genes used to estimate the biological correlation} 
 }
 \seealso{
\code{\link{lmFit}}, \code{\link{eBayes}}, \code{\link{contrasts.fit}}
}

\author{Belinda Phipson and Gordon Smyth}

\references{
Majewski, IJ, Ritchie, ME, Phipson, B, Corbin, J, Pakusch, M, Ebert, A, Busslinger, M, Koseki, H, Hu, Y, Smyth, GK, Alexander, WS, Hilton, DJ, and Blewitt, ME (2010).
Opposing roles of polycomb repressive complexes in hematopoietic stem and progenitor cells.
\emph{Blood} 116, 731-739.
\url{http://www.bloodjournal.org/content/116/5/731}

Phipson, B. (2013).
\emph{Empirical Bayes modelling of expression profiles and their associations}.
PhD Thesis. University of Melbourne, Australia.
\url{http://repository.unimelb.edu.au/10187/17614}

Ritchie, ME, Phipson, B, Wu, D, Hu, Y, Law, CW, Shi, W, and Smyth, GK (2015).
limma powers differential expression analyses for RNA-sequencing and microarray studies.
\emph{Nucleic Acids Research} 43, e47.
\url{http://nar.oxfordjournals.org/content/43/7/e47}
}

\examples{
#  Simulate gene expression data

#  Three conditions (Control, A and B) and 1000 genes
ngene <- 1000
mu.A <- mu.B <- mu.ctrl <- rep(5,ngene)

#  200 genes are differentially expressed.
#  All are up in condition A and down in B
#  so the biological correlation is negative.
mu.A[1:200] <- mu.ctrl[1:200]+2
mu.B[1:200] <- mu.ctrl[1:200]-2

#  Two microarrays for each condition
mu <- cbind(mu.ctrl,mu.ctrl,mu.A,mu.A,mu.B,mu.B)
y <- matrix(rnorm(6000,mean=mu,sd=1),ngene,6)

# two experimental groups and one control group with two replicates each
group <- factor(c("Ctrl","Ctrl","A","A","B","B"), levels=c("Ctrl","A","B"))
design <- model.matrix(~group)

# fit a linear model
fit <- lmFit(y,design)
fit <- eBayes(fit)

#  Estimate biological correlation between the logFC profiles
#  for A-vs-Ctrl and B-vs-Ctrl
genas(fit, coef=c(2,3), plot=TRUE, subset="F")
}
