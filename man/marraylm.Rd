\name{MArrayLM-class}
\docType{class}
\alias{MArrayLM-class}
\title{Microarray Linear Model Fit - class}

\description{
A list-based S4 class for storing the results of fitting gene-wise linear models to a set of microarrays.
Objects are normally created by \code{\link{lmFit}}, and additional components are added by \code{\link{eBayes}}.
}

\section{Components}{
  \code{MArrayLM} objects do not contain any slots (apart from \code{.Data}) but they should contain the following list components:
  \tabular{ll}{
    {\code{coefficients}} \tab {matrix containing fitted coefficients or contrasts}\cr
    {\code{stdev.unscaled}} \tab {matrix containing unscaled standard deviations of the coefficients or contrasts}\cr
    {\code{sigma}} \tab {numeric vector containing residual standard deviations for each gene}\cr
    {\code{df.residual}} \tab {numeric vector containing residual degrees of freedom for each gene}
  }
  The following additional components may be created by \code{lmFit}:
  \tabular{ll}{
    {\code{Amean}} \tab {numeric vector containing the average log-intensity for each probe over all the arrays in the original linear model fit. Note this vector does not change when a contrast is applied to the fit using \code{contrasts.fit}.}\cr
    {\code{genes}} \tab {data.frame containing probe annotation.}\cr
    {\code{design}} \tab {design matrix.}\cr
    {\code{cov.coefficients}} \tab {numeric matrix giving the unscaled covariance matrix of the estimable coefficients}\cr
    {\code{pivot}} \tab {integer vector giving the order of coefficients in \code{cov.coefficients}. Is computed by the QR-decomposition of the design matrix.}\cr
    \code{qr} \tab QR-decomposition of the design matrix (if the fit involved no weights or missing values).\cr
    \dots \tab other components returned by \code{lm.fit} (if the fit involved no weights or missing values).
  }
  The following component may be added by \code{contrasts.fit}:
  \tabular{ll}{
    {\code{contrasts}} \tab {numeric matrix defining contrasts of coefficients for which results are desired.}
  }
  The following components may be added by \code{eBayes}:
  \tabular{ll}{
    {\code{s2.prior}} \tab {numeric value giving empirical Bayes estimated prior value for residual variances}\cr
    {\code{df.prior}} \tab {numeric vector giving empirical Bayes estimated degrees of freedom associated with \code{s2.prior} for each gene}\cr
    {\code{df.total}} \tab {numeric vector giving total degrees of freedom used for each gene, usually equal to \code{df.prior + df.residual}.}\cr
    {\code{s2.post}} \tab {numeric vector giving posterior residual variances}\cr
    {\code{var.prior}} \tab {numeric vector giving empirical Bayes estimated prior variance for each true coefficient}\cr
    {\code{F}} \tab {numeric vector giving moderated F-statistics for testing all contrasts equal to zero}\cr
    {\code{F.p.value}} \tab {numeric vector giving p-value corresponding to \code{F.stat}}\cr
    {\code{t}} \tab {numeric matrix containing empirical Bayes t-statistics}\cr
  }
}

\section{Methods}{
\code{MArrayLM} objects will return dimensions and hence functions such as \code{\link[limma:dim]{dim}}, \code{\link[base:nrow]{nrow}} and \code{\link[base:nrow]{ncol}} are defined. 
\code{MArrayLM} objects inherit a \code{show} method from the virtual class \code{LargeDataObject}.

The functions \code{\link{eBayes}}, \code{\link{decideTests}} and \code{\link{classifyTestsF}} accept \code{MArrayLM} objects as arguments.
}

\author{Gordon Smyth}

\seealso{
  \link{02.Classes} gives an overview of all the classes defined by this package.
}

\keyword{classes}
\keyword{regression}
