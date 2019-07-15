\title{Plot of regularized linear discriminant functions for microarray data}
\name{plotRLDF}
\alias{plotRLDF}
\description{
Plot regularized linear discriminant functions for classifying samples based on expression data.
}
\usage{
plotRLDF(y, design = NULL, z = NULL, nprobes = 100, plot = TRUE,
         labels.y = NULL, labels.z = NULL, pch.y = NULL, pch.z = NULL,
         col.y = "black", col.z = "black",
         show.dimensions = c(1,2), ndim = max(show.dimensions),
         var.prior = NULL, df.prior = NULL, trend = FALSE, robust = FALSE, \dots)
}
\arguments{
 \item{y}{the training dataset. Can be any data object which can be coerced to a matrix, such as \code{ExpressionSet} or \code{EList}.}
 \item{design}{design matrix defining the training groups to be distinguished. The first column is assumed to represent the intercept.
 Defaults to \code{model.matrix(~factor(labels.y))}.}
 \item{z}{the dataset to be classified.  Can be any data object which can be coerced to a matrix, such as \code{ExpressionSet} or \code{EList}.  Rows must correspond to rows of \code{y}.}
 \item{nprobes}{number of probes to be used for the calculations. The probes will be selected by moderated F statistic.}
 \item{plot}{logical, should a plot be created?}
 \item{labels.y}{character vector of sample names or labels in \code{y}. Defaults to \code{colnames(y)} or failing that to \code{1:n}.}
 \item{labels.z}{character vector of sample names or labels in \code{z}. Defaults to \code{colnames(z)} or failing that to \code{letters[1:n]}.}
 \item{pch.y}{plotting symbol or symbols for \code{y}. See \code{\link{points}} for possible values. Takes precedence over \code{labels.y} if both are specified.}
 \item{pch.z}{plotting symbol or symbols for \code{y}. See \code{\link{points}} for possible values. Takes precedence over \code{labels.z} if both are specified.}
 \item{col.y}{colors for the plotting \code{labels.y}.}
 \item{col.z}{colors for the plotting \code{labels.z}.}
 \item{show.dimensions}{integer vector of length two indicating which two discriminant functions to plot. Functions are in decreasing order of discriminatory power.}
 \item{ndim}{number of discriminant functions to compute}
 \item{var.prior}{prior variances, for regularizing the within-group covariance matrix. By default is estimated by \code{squeezeVar}.}
 \item{df.prior}{prior degrees of freedom for regularizing the within-group covariance matrix. By default is estimated by \code{squeezeVar}.}
 \item{trend}{logical, should a trend be estimated for \code{var.prior}?  See \code{eBayes} for details.  Only used if \code{var.prior} or \code{df.prior} are \code{NULL}.}
 \item{robust}{logical, should \code{var.prior} and \code{df.prior} be estimated robustly?  See \code{eBayes} for details.  Only used if \code{var.prior} or \code{df.prior} are \code{NULL}.}
 \item{\dots}{any other arguments are passed to \code{plot}.}
}

\details{
The function builds discriminant functions from the training data (\code{y}) and applies them to the test data (\code{z}).
The method is a variation on classifical linear discriminant functions (LDFs), in that the within-group covariance matrix is regularized to ensure that it is invertible, with eigenvalues bounded away from zero.
The within-group covariance matrix is squeezed towards a diagonal matrix with empirical Bayes posterior variances as diagonal elements.

The calculations are based on a filtered list of probes.
The \code{nprobes} probes with largest moderated F statistics are used to discriminate.

The \code{ndim} argument allows all required LDFs to be computed even though only two are plotted.
}

\note{
The default values for \code{df.prior} and \code{var.prior} were changed in limma 3.27.10.
Previously these were preset values.
Now the default is to estimate them using \code{squeezeVar}.
}

\value{
If \code{plot=TRUE} a plot is created on the current graphics device.
A list containing the following components is (invisibly) returned:
  \item{training}{numeric matrix with \code{ncol(y)} rows and \code{ndim} columns containing discriminant functions evaluated for the training data.}
  \item{predicting}{numeric matrix with \code{ncol(z)} rows and \code{ndim} columns containing discriminant functions evalulated on the classification data.}
  \item{top}{integer vector of length \code{nprobes} giving indices of probes used.}
  \item{metagenes}{numeric matrix with \code{nprobes} rows and \code{ndim} columns containing probe weights defining each discriminant function.}
  \item{singular.values}{singular.values showing the predictive power of each discriminant function.}
  \item{rank}{maximum number of discriminant functions with singular.values greater than zero.}
  \item{var.prior}{numeric vector of prior variances.}
  \item{df.prior}{numeric vector of prior degrees of freedom.}
}

\author{Gordon Smyth, Di Wu and Yifang Hu}

\seealso{
\code{lda} in package \code{MASS}
}

\examples{
# Simulate gene expression data for 1000 probes and 6 microarrays.
# Samples are in two groups
# First 50 probes are differentially expressed in second group
sd <- 0.3*sqrt(4/rchisq(1000,df=4))
y <- matrix(rnorm(1000*6,sd=sd),1000,6)
rownames(y) <- paste("Gene",1:1000)
y[1:50,4:6] <- y[1:50,4:6] + 2

z <- matrix(rnorm(1000*6,sd=sd),1000,6)
rownames(z) <- paste("Gene",1:1000)
z[1:50,4:6] <- z[1:50,4:6] + 1.8
z[1:50,1:3] <- z[1:50,1:3] - 0.2

design <- cbind(Grp1=1,Grp2vs1=c(0,0,0,1,1,1))
options(digit=3)

# Samples 1-6 are training set, samples a-f are test set:
plotRLDF(y, design, z=z, col.y="black", col.z="red")
legend("top", pch=16, col=c("black","red"), legend=c("Training","Predicted"))
}

\keyword{hplot}

