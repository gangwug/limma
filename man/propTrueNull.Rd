\name{propTrueNull}
\alias{convest}
\alias{propTrueNull}
\alias{estimate.m0}

\title{Estimate Proportion of True Null Hypotheses}

\description{
Estimate the proportion of true null hypotheses from a vector of p-values.
}

\usage{
propTrueNull(p, method="lfdr", nbins=20, \dots)
convest(p, niter=100, plot=FALSE, report=FALSE, file="", tol=1e-6)
}

\arguments{
  \item{p}{numeric vector of p-values.}
  \item{method}{estimation method. Choices are \code{"lfdr"}, \code{"mean"}, \code{"hist"} or \code{"convest"}.}
  \item{nbins}{number of histogram bins (if \code{method="hist"}).}
  \item{niter}{number of iterations to be used in fitting the convex, decreasing density for the p-values.}
  \item{plot}{logical, should updated plots of fitted convex decreasing p-value density be produced at each iteration?}
  \item{report}{logical, should the estimated proportion be printed at each iteration?}
  \item{file}{name of file to which to write the report. Defaults to standard output.}
  \item{tol}{accuracy of the bisectional search for finding a new convex combination of the current iterate and the mixing density}
  \item{\dots}{other arguments are passed to \code{convest} if \code{method="convest"}.}
}

\details{
The proportion of true null hypotheses in a collection of hypothesis tests is often denoted pi0.
This function estimates pi0 from a vector of p-values.

\code{method="lfdr"} implements the method of Phipson (2013) based on averaging local false discovery rates across the p-values.

\code{method="mean"} is a very simple method based on averaging the p-values.  It gives a slightly smaller estimate than \code{2*mean(p)}.

\code{method="hist"} implements the histogram method of Mosig et al (2001) and Nettleton et al (2006).

\code{method="convest"} calls \code{convest}, which implements the method of Langaas et al (2005) based on a convex decreasing density estimate.
}

\value{
Numeric value in the interval [0,1] representing the estimated proportion of true null hypotheses.
}

\references{
Langaas, M, Ferkingstad, E, and Lindqvist, B (2005).
Estimating the proportion of true null hypotheses, with application to DNA microarray data.
\emph{Journal of the Royal Statistical Society Series} B 67, 555-572.

Mosig MO, Lipkin E, Khutoreskaya G, Tchourzyna E, Soller M, Friedmann A (2001).
A whole genome scan for quantitative trait loci affecting milk protein percentage in Israeli-Holstein cattle, by means of selective milk DNA pooling in a daughter design, using an adjusted false discovery rate criterion.
\emph{Genetics} 157, 1683-1698.

Nettleton D, Hwang JTG, Caldo RA, Wise RP (2006).
Estimating the number of true null hypotheses from a histogram of p values.
\emph{Journal of Agricultural, Biological, and Environmental Statistics} 11, 337-356.

Phipson, B (2013).
Empirical Bayes Modelling of Expression Profiles and Their Associations.
PhD Thesis, University of Melbourne, Australia.
\url{http://repository.unimelb.edu.au/10187/17614}

Ritchie, ME, Phipson, B, Wu, D, Hu, Y, Law, CW, Shi, W, and Smyth, GK (2015).
limma powers differential expression analyses for RNA-sequencing and microarray studies.
\emph{Nucleic Acids Research} 43, e47.
\url{http://nar.oxfordjournals.org/content/43/7/e47}
}

\author{Belinda Phipson and Gordon Smyth for \code{propTrueNull}. Egil Ferkingstad, Mette Langaas and Marcus Davy for \code{convest}.}

\seealso{
See \link{08.Tests} for other functions for producing or interpreting p-values.
}

\examples{
# Test statistics
z <- rnorm(200)

# First 40 are have non-zero means
z[1:40] <- z[1:40]+2

# True pi0
160/200

# Two-sided p-values
p <- 2*pnorm(-abs(z))

# Estimate pi0
propTrueNull(p, method="lfdr")
propTrueNull(p, method="hist")
}

\keyword{htest}
