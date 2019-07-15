\name{lmFit}
\alias{lmFit}
\title{Linear Model for Series of Arrays}
\description{Fit linear model for each gene given a series of arrays}
\usage{
lmFit(object, design=NULL, ndups=1, spacing=1, block=NULL, correlation, weights=NULL,
      method="ls", \dots)
}
\arguments{
  \item{object}{A matrix-like data object containing log-ratios or log-expression values for a series of arrays, with rows corresponding to genes and columns to samples.
  Any type of data object that can be processed by \code{\link{getEAWP}} is acceptable.}
  \item{design}{the design matrix of the microarray experiment, with rows corresponding to arrays and columns to coefficients to be estimated.  Defaults to the unit vector meaning that the arrays are treated as replicates.} 
  \item{ndups}{positive integer giving the number of times each distinct probe is printed on each array.}
  \item{spacing}{positive integer giving the spacing between duplicate occurrences of the same probe, \code{spacing=1} for consecutive rows.}
  \item{block}{vector or factor specifying a blocking variable on the arrays. Has length equal to the number of arrays. Must be \code{NULL} if \code{ndups>2}.}
  \item{correlation}{the inter-duplicate or inter-technical replicate correlation}
  \item{weights}{non-negative precision weights.  Can be a numeric matrix of individual weights of same size as the object expression matrix, or a numeric vector of array weights with length equal to \code{ncol} of the expression matrix, or a numeric vector of gene weights with length equal to \code{nrow} of the expression matrix.}
  \item{method}{fitting method; \code{"ls"} for least squares or \code{"robust"} for robust regression}
  \item{\dots}{other optional arguments to be passed to \code{lm.series}, \code{gls.series} or \code{mrlm}}
}

\value{
An \code{\link[limma:marraylm]{MArrayLM}} object containing the result of the fits.

The rownames of \code{object} are preserved in the fit object and can be retrieved by \code{rownames(fit)} where \code{fit} is output from \code{lmFit}.
The column names of \code{design} are preserved as column names and can be retrieved by \code{colnames(fit)}.
}

\details{
This function fits multiple linear models by weighted or generalized least squares.
It accepts data from a experiment involving a series of microarrays with the same set of probes.
A linear model is fitted to the expression data for each probe.
The expression data should be log-ratios for two-color array platforms or log-expression values for one-channel platforms.
(To fit linear models to the individual channels of two-color array data, see \code{\link{lmscFit}}.)
The coefficients of the fitted models describe the differences between the RNA sources hybridized to the arrays.
The probe-wise fitted model results are stored in a compact form suitable for further processing by other functions in the limma package.

The function allows for missing values and accepts quantitative precision weights through the \code{weights} argument.
It also supports two different correlation structures.
If \code{block} is not \code{NULL} then different arrays are assumed to be correlated.
If \code{block} is \code{NULL} and \code{ndups} is greater than one then replicate spots on the same array are assumed to be correlated.  
It is not possible at this time to fit models with both a block structure and a duplicate-spot correlation structure simultaneously.

If \code{object} is a matrix then it should contain log-ratios or log-expression data with rows corresponding to probes and columns to arrays.
(A numeric vector is treated the same as a matrix with one column.)
For objects of other classes, a matrix of expression values is taken from the appropriate component or slot of the object.
If \code{object} is of class \code{MAList} or \code{marrayNorm}, then the matrix of log-ratios (M-values) is extracted.
If \code{object} is of class \code{ExpressionSet}, then the expression matrix is extracted.
(This may contain log-expression or log-ratio values, depending on the platform.)
If \code{object} is of class \code{PLMset} then the matrix of chip coefficients \code{chip.coefs} is extracted.

The arguments \code{design}, \code{ndups}, \code{spacing} and \code{weights} will be extracted from the data \code{object} if available and do not normally need to set explicitly in the call.
On the other hand, if any of these are set in the function call then they will over-ride the slots or components in the data \code{object}.
If \code{object} is an \code{PLMset}, then weights are computed as \code{1/pmax(object@se.chip.coefs, 1e-05)^2}.
If \code{object} is an \code{ExpressionSet} object, then weights are not computed.

If the argument \code{block} is used, then it is assumed that \code{ndups=1}.

The \code{correlation} argument has a default value of \code{0.75}, but in normal use this default value should not be relied on and the correlation value should be estimated using the function \code{duplicateCorrelation}.
The default value is likely to be too high in particular if used with the \code{block} argument.

The actual linear model computations are done by passing the data to one the lower-level functions \code{lm.series}, \code{gls.series} or \code{mrlm}.
The function \code{mrlm} is used if \code{method="robust"}.
If \code{method="ls"}, then \code{gls.series} is used if a correlation structure has been specified, i.e., if \code{ndups>1} or \code{block} is non-null and \code{correlation} is different from zero.
If \code{method="ls"} and there is no correlation structure, \code{lm.series} is used.
}

\seealso{
\code{lmFit} uses \code{\link{getEAWP}} to extract expression values, gene annotation and so from the data \code{object}.

An overview of linear model functions in limma is given by \link{06.LinearModels}.
}

\author{Gordon Smyth}
\examples{
# Simulate gene expression data for 100 probes and 6 microarrays
# Microarray are in two groups
# First two probes are differentially expressed in second group
# Std deviations vary between genes with prior df=4
sd <- 0.3*sqrt(4/rchisq(100,df=4))
y <- matrix(rnorm(100*6,sd=sd),100,6)
rownames(y) <- paste("Gene",1:100)
y[1:2,4:6] <- y[1:2,4:6] + 2
design <- cbind(Grp1=1,Grp2vs1=c(0,0,0,1,1,1))
options(digits=3)

# Ordinary fit
fit <- lmFit(y,design)
fit <- eBayes(fit)
topTable(fit,coef=2)
dim(fit)
colnames(fit)
rownames(fit)[1:10]
names(fit)

# Fold-change thresholding
fit2 <- treat(fit,lfc=0.1)
topTreat(fit2,coef=2)

# Volcano plot
volcanoplot(fit,coef=2,highlight=2)

# Mean-difference plot
plotMD(fit,column=2)

# Q-Q plot of moderated t-statistics
qqt(fit$t[,2],df=fit$df.residual+fit$df.prior)
abline(0,1)

# Various ways of writing results to file
\dontrun{write.fit(fit,file="exampleresults.txt")}
\dontrun{write.table(fit,file="exampleresults2.txt")}

# Fit with correlated arrays
# Suppose each pair of arrays is a block
block <- c(1,1,2,2,3,3)
dupcor <- duplicateCorrelation(y,design,block=block)
dupcor$consensus.correlation
fit3 <- lmFit(y,design,block=block,correlation=dupcor$consensus)

# Fit with duplicate probes
# Suppose two side-by-side duplicates of each gene
rownames(y) <- paste("Gene",rep(1:50,each=2))
dupcor <- duplicateCorrelation(y,design,ndups=2)
dupcor$consensus.correlation
fit4 <- lmFit(y,design,ndups=2,correlation=dupcor$consensus)
dim(fit4)
fit4 <- eBayes(fit4)
topTable(fit4,coef=2)
}
\keyword{models}
\keyword{regression}
