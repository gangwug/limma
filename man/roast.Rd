\name{roast}
\alias{roast}
\alias{roast.default}
\alias{mroast}
\alias{mroast.default}
\alias{Roast-class}
\alias{show,Roast-method}
\alias{fry}
\alias{fry.default}
\title{Rotation Gene Set Tests}
\description{
Rotation gene set testing for linear models.
}

\usage{
\method{roast}{default}(y, index = NULL, design = NULL, contrast = ncol(design), geneid = NULL,
      set.statistic = "mean", gene.weights = NULL, var.prior = NULL, df.prior = NULL,
      nrot = 999, approx.zscore = TRUE, \dots)
\method{mroast}{default}(y, index = NULL, design = NULL, contrast = ncol(design), geneid = NULL,
       set.statistic = "mean", gene.weights = NULL, var.prior = NULL, df.prior = NULL,
       nrot = 999, approx.zscore = TRUE, adjust.method = "BH",
       midp = TRUE, sort = "directional", \dots)
\method{fry}{default}(y, index = NULL, design = NULL, contrast = ncol(design), geneid = NULL,
      standardize = "posterior.sd", sort = "directional", \dots)
}

\arguments{
  \item{y}{numeric matrix giving log-expression or log-ratio values for a series of microarrays, or any object that can coerced to a matrix including \code{ExpressionSet}, \code{MAList}, \code{EList} or \code{PLMSet} objects.
        Rows correspond to probes and columns to samples.
        If either \code{var.prior} or \code{df.prior} are \code{NULL}, then \code{y} should contain values for all genes on the arrays. If both prior parameters are given, then only \code{y} values for the test set are required.}
  \item{index}{index vector specifying which rows (probes) of \code{y} are in the test set.
        Can be a vector of integer indices, or a logical vector of length \code{nrow(y)}, or a vector of gene IDs corresponding to entries in \code{geneid}.
        Alternatively it can be a data.frame with the first column containing the index vector and the second column containing directional gene contribution weights.
        For \code{mroast} or \code{fry}, \code{index} is a list of index vectors or a list of data.frames. }
  \item{design}{design matrix}
  \item{contrast}{contrast for which the test is required.
        Can be an integer specifying a column of \code{design}, or the name of a column of \code{design}, or a numeric contrast vector of length equal to the number of columns of \code{design}.}
  \item{geneid}{gene identifiers corresponding to the rows of \code{y}.
        Can be either a vector of length \code{nrow(y)} or the name of the column of \code{y$genes} containing the gene identifiers.
        Defaults to \code{rownames(y)}.}
  \item{set.statistic}{summary set statistic. Possibilities are \code{"mean"},\code{"floormean"},\code{"mean50"} or \code{"msq"}.}
  \item{gene.weights}{numeric vector of directional (positive or negative) contribution weights specifying the size and direction of the contribution of each probe to the gene set statistics.
        For \code{mroast}, this vector must have length equal to \code{nrow(y)}.
        For \code{roast}, can be of length \code{nrow(y)} or of length equal to the number of genes in the test set.} 
  \item{var.prior}{prior value for residual variances. If not provided, this is estimated from all the data using \code{squeezeVar}.}
  \item{df.prior}{prior degrees of freedom for residual variances. If not provided, this is estimated using \code{squeezeVar}.}
  \item{nrot}{number of rotations used to compute the p-values.}
  \item{approx.zscore}{logical, if \code{TRUE} then a fast approximation is used to convert t-statistics into z-scores prior to computing set statistics. If \code{FALSE}, z-scores will be exact.}
  \item{adjust.method}{method used to adjust the p-values for multiple testing. See \code{\link{p.adjust}} for possible values.}
  \item{midp}{logical, should mid-p-values be used in instead of ordinary p-values when adjusting for multiple testing?}
  \item{sort}{character, whether to sort output table by directional p-value (\code{"directional"}), non-directional p-value (\code{"mixed"}), or not at all (\code{"none"}).}
  \item{standardize}{how to standardize for unequal probewise variances. Possibilities are \code{"residual.sd"}, \code{"posterior.sd"} or \code{"none"}.}
  \item{\dots}{any argument that would be suitable for \code{\link{lmFit}} or \code{\link{eBayes}} can be included.}
}

\value{
\code{roast} produces an object of class \code{"Roast"}.
This consists of a list with the following components:
  \item{p.value}{data.frame with columns \code{Active.Prop} and \code{P.Value}, giving the proportion of genes in the set contributing materially to significance and estimated p-values, respectively.
Rows correspond to the alternative hypotheses Down, Up, UpOrDown (two-sided) and Mixed.}
  \item{var.prior}{prior value for residual variances.}
  \item{df.prior}{prior degrees of freedom for residual variances.}

\code{mroast} produces a data.frame with a row for each set and the following columns:
  \item{NGenes}{number of genes in set}
	\item{PropDown}{proportion of genes in set with \code{z < -sqrt(2)}}
	\item{PropUp}{proportion of genes in set with \code{z > sqrt(2)}}
	\item{Direction}{direction of change, \code{"Up"} or \code{"Down"}}
	\item{PValue}{two-sided directional p-value}
	\item{FDR}{two-sided directional false discovery rate}
	\item{PValue.Mixed}{non-directional p-value}
	\item{FDR.Mixed}{non-directional false discovery rate}

\code{fry} produces the same output format as \code{mroast} but without the columns \code{PropDown} and \code{ProbUp}.
}

\details{
These functions implement the ROAST gene set tests proposed by Wu et al (2010).
They perform \emph{self-contained} gene set tests in the sense defined by Goeman and Buhlmann (2007).
For \emph{competitive} gene set tests, see \code{\link{camera}}.
For a gene set enrichment analysis style analysis using a database of gene sets, see \code{\link{romer}}.

\code{roast} and \code{mroast} test whether any of the genes in the set are differentially expressed.
They can be used for any microarray experiment that can be represented by a linear model.
The design matrix for the experiment is specified as for the \code{\link{lmFit}} function, and the contrast of interest is specified as for the \code{\link{contrasts.fit}} function.
This allows users to focus on differential expression for any coefficient or contrast in a linear model.
If \code{contrast} is not specified, then the last coefficient in the linear model will be tested.

The argument \code{index} is often made using \link{ids2indices} but does not have to be.
Each set to be tested is represented by a vector of row numbers or a vector of gene IDs.
Gene IDs should correspond to either the rownames of \code{y} or the entries of \code{geneid}.

The argument \code{gene.weights} allows directional contribution weights to be set for individual genes in the set.
This is often useful, because it allows each gene to be flagged as to its direction and magnitude of change based on prior experimentation.
A typical use is to make the \code{gene.weights} \code{1} or \code{-1} depending on whether the gene is up or down-regulated in the pathway under consideration.
Probes with directional weights of opposite signs are expected to have expression changes in opposite directions.
If there are multiple sets to be tested, then set-specific gene weights can be included as part of the \code{index}.
If any of the entries of \code{index} are data.frames, then the second column will be assumed to be gene contribution weights.
All three functions (\code{roast}, \code{mroast} and \code{fry}) support set-specific gene contribution weights as part of an \code{index} data.frame.

Note that the contribution weights set by \code{gene.weights} are different in nature and purpose to the precision weights set by the \code{weights} argument to \code{lmFit}.
\code{gene.weights} control the contribution of each gene to the formation of the gene set statistics, and can be positive or negative.
\code{weights} indicate the precision of the expression measurements and should be positive.
The \code{weights} are used to construct genewise test statistics whereas \code{gene.weights} are used to combine the genewise test statistics.

The arguments \code{df.prior} and \code{var.prior} have the same meaning as in the output of the \code{\link{eBayes}} function.
If these arguments are not supplied, then they are estimated exactly as is done by \code{eBayes}.

The gene set statistics \code{"mean"}, \code{"floormean"}, \code{"mean50"} and \code{msq} are defined by Wu et al (2010).
The different gene set statistics have different sensitivities to small number of genes.
If \code{set.statistic="mean"} then the set will be statistically significantly only when the majority of the genes are differentially expressed.
\code{"floormean"} and \code{"mean50"} will detect as few as 25\% differentially expressed.
\code{"msq"} is sensitive to even smaller proportions of differentially expressed genes, if the effects are reasonably large.

The output gives p-values three possible alternative hypotheses, 
\code{"Up"} to test whether the genes in the set tend to be up-regulated, with positive t-statistics,
\code{"Down"} to test whether the genes in the set tend to be down-regulated, with negative t-statistics,
and \code{"Mixed"} to test whether the genes in the set tend to be differentially expressed, without regard for direction.

\code{roast} estimates p-values by simulation, specifically by random rotations of the orthogonalized residuals (Langsrud, 2005), so p-values will vary slightly from run to run.
To get more precise p-values, increase the number of rotations \code{nrot}.
The p-value is computed as \code{(b+1)/(nrot+1)} where \code{b} is the number of rotations giving a more extreme statistic than that observed (Phipson and Smyth, 2010).
This means that the smallest possible p-value is \code{1/(nrot+1)}.

\code{mroast} does roast tests for multiple sets, including adjustment for multiple testing.
By default, \code{mroast} reports ordinary p-values but uses mid-p-values (Routledge, 1994) at the multiple testing stage.
Mid-p-values are probably a good choice when using false discovery rates (\code{adjust.method="BH"}) but not when controlling the family-wise type I error rate (\code{adjust.method="holm"}).

\code{fry} is a fast approximation to \code{mroast}.
In the special case that \code{df.prior} is large and \code{set.statistic="mean"}, \code{fry} gives the same result as \code{mroast} with an infinite number of rotations.
In other circumstances, when genes have different variances, \code{fry} uses a standardization strategy to approximate the \code{mroast} results.
Using \code{fry} may be advisable when performing tests for a large number of sets, because it is fast and because the \code{fry} p-values are not limited by the number of rotations performed.
}

\note{
The default setting for the set statistic was changed in limma 3.5.9 (3 June 2010) from \code{"msq"} to \code{"mean"}.
}

\seealso{
See \link{10.GeneSetTests} for a description of other functions used for gene set testing.
}
\author{Gordon Smyth and Di Wu}

\references{
Goeman, JJ, and Buhlmann, P (2007).
Analyzing gene expression data in terms of gene sets: methodological issues.
\emph{Bioinformatics} 23, 980-987. 

Langsrud, O (2005).
Rotation tests.
\emph{Statistics and Computing} 15, 53-60.

Phipson B, and Smyth GK (2010).
Permutation P-values should never be zero: calculating exact P-values when permutations are randomly drawn.
\emph{Statistical Applications in Genetics and Molecular Biology}, Volume 9, Article 39.
\url{http://www.statsci.org/smyth/pubs/PermPValuesPreprint.pdf}

Routledge, RD (1994).
Practicing safe statistics with the mid-p.
\emph{Canadian Journal of Statistics} 22, 103-110.

Wu, D, Lim, E, Francois Vaillant, F, Asselin-Labat, M-L, Visvader, JE, and Smyth, GK (2010). ROAST: rotation gene set tests for complex microarray experiments.
\emph{Bioinformatics} 26, 2176-2182.
\url{http://bioinformatics.oxfordjournals.org/content/26/17/2176}
}

\examples{
y <- matrix(rnorm(100*4),100,4)
design <- cbind(Intercept=1,Group=c(0,0,1,1))

# First set of 5 genes contains 3 that are genuinely differentially expressed
index1 <- 1:5
y[index1,3:4] <- y[index1,3:4]+3

# Second set of 5 genes contains none that are DE
index2 <- 6:10

roast(y,index1,design,contrast=2)
fry(y,list(set1=index1,set2=index2),design,contrast=2)
}

\keyword{gene set test}
\concept{gene set test}
\concept{gene set enrichment analysis}
