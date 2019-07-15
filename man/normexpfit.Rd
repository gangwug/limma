\name{normexp.fit}
\alias{normexp.fit}
\title{Fit Normal+Exp Convolution Model to Observed Intensities}
\description{
Fit the normal+exponential convolution model to a vector of observed intensities.
The normal part represents the background and the exponential part represents the signal intensities.
This function is called by \code{backgroundCorrect} and is not normally called directly by users.
}
\usage{
normexp.fit(x, method="saddle", n.pts=NULL, trace=FALSE)
}
\arguments{
  \item{x}{numeric vector of (background corrected) intensities}
  \item{method}{method used to estimate the three parameters. Choices for \code{normexp.fit} are \code{"mle"}, \code{"saddle"}, \code{"rma"} and \code{"rma75"}.}
  \item{n.pts}{number of quantiles of \code{x} to use for the fit. If \code{NULL} then all values of \code{x} will be used.}
  \item{trace}{logical, if \code{TRUE}, tracing information on the progress of the optimization is given.}
}
\details{
The Normal+Exp (normexp) convolution model is a mathematical model representing microarray intensity data for the purposes of background correction.
It was proposed originally as part of the RMA algorithm for Affymetrix microarray data.
For two-color microarry data, the normexp background correction method was introduced and compared with other methods by Ritchie et al (2007).

This function uses maximum likelihood estimation to fit the normexp model to background-corrected intensities.
The model assumes that the observed intensities are the sum of background and signal components,
the background being normal and the signal being exponential distributed.

The likelihood may be computed exactly (\code{method="mle"}) or approximated using a saddle-point approximation (\code{method="saddle"}).
The saddle-point approximation was proposed by Ritchie et al (2007).
Silver et al (2008) added some computational refinements to the saddle-point approximation, making it more reliable in practice, and developed the exact likelihood maximization algorithm.
The \code{"mle"} method uses the best performing algorithm from Silver et al (2008), which
calls the optimization function \code{\link{nlminb}} with analytic first and second derivatives.
Derivatives are computed with respect to the normal-mean, the log-normal-variance and the log-exponential-mean.

Two ad-hoc estimators are also available which do not require iterative estimation.
\code{"rma"} results in a call to the \code{bg.parameters} function of the affy package.
This provides the kernel estimation method that is part of the RMA algorithm for Affymetrix data.
\code{"rma75"} uses the similar but less biased RMA-75 method from McGee and Chen (2006).

If the length \code{x} is very large, it may be worth saving computation time by setting \code{n.pts} to a value less than the total number of probes, for example \code{n.pts=2^14}.
}
\value{
A list containing the components
  \item{par}{numeric vector giving estimated values of the mean and log-standard-deviation of the background-normal part and the log-mean of the signal-exponential part.}
  \item{m2loglik}{numeric scalar giving minus twice the maximized log-likelihood}
  \item{convergence}{integer code indicating successful convergence or otherwise of the optimization.}
}
\author{Gordon Smyth and Jeremy Silver}

\references{
McGee, M., and Chen, Z. (2006).
Parameter estimation for the exponential-normal convolution model for background correction of Affymetrix GeneChip data.
\emph{Stat Appl Genet Mol Biol}, 5(1), Article 24.

Ritchie, M. E., Silver, J., Oshlack, A., Silver, J., Holmes, M., Diyagama, D., Holloway, A., and Smyth, G. K. (2007).
A comparison of background correction methods for two-colour microarrays.
\emph{Bioinformatics} 23, 2700-2707.
\url{http://bioinformatics.oxfordjournals.org/content/23/20/2700}

Silver, JD, Ritchie, ME, and Smyth, GK (2009).
Microarray background correction: maximum likelihood estimation for the normal-exponential convolution.
\emph{Biostatistics} 10, 352-363.
\url{http://biostatistics.oxfordjournals.org/content/10/2/352}
}

\seealso{
\code{\link{normexp.signal}}, \code{\link{normexp.fit.control}}.
Also \code{\link[affy:bg.adjust]{bg.parameters} in the affy package.
}

An overview of background correction functions is given in \code{\link{04.Background}}.
}
\examples{
x <- c(2,3,1,10,3,20,5,6)
out <- normexp.fit(x)
normexp.signal(out$par, x=x)
}
\keyword{models}
