\name{weightedLoess}
\alias{weightedLowess}

\title{Lowess fit with weighting}

\description{Fit robust lowess curves of degree 1 to weighted covariates and responses.}

\usage{
weightedLowess(x, y, weights = rep(1, length(y)),
               delta=NULL, npts = 200, span = 0.3, iterations = 4)
}

\arguments{
\item{x}{a numeric vector of covariates}
\item{y}{a numeric vector of response values}
\item{weights}{a numeric vector containing frequency weights for each covariate}
\item{delta}{a numeric scalar specifying the maximum distance between adjacent points}
\item{npts}{an integer scalar specifying the approximate number of points to use when computing \code{delta}}
\item{span}{a numeric scalar specifying the width of the smoothing window as a proportion of the total weight}
\item{iterations}{an integer scalar specifying the number of robustifying iterations}
}

\value{
A list of numeric vectors for the fitted responses, the residuals, the robustifying weights and the chosen delta.
}

\details{
This function extends the lowess algorithm to handle non-negative prior weights. These weights are
used during span calculations such that the span distance for each point must include the
specified proportion of all weights. They are also applied during weighted linear regression to 
compute the fitted value (in addition to the tricube weights determined by \code{span}). For integer 
weights, the prior weights are equivalent to using \code{rep(\dots, w)} on \code{x} and \code{y} prior to fitting.

For large vectors, running time is reduced by only performing locally weighted regression for 
several points. Fitted values for all points adjacent to the chosen points are computed by 
linear interpolation between the chosen points. For this purpose, the first and last points are always
chosen. Note that the regression itself uses all (neighbouring) points.

Points are defined as adjacent to a chosen point if the distance to the latter is positive and less 
than \code{delta}. The first chosen point is that corresponding to the smallest covariate; the
next chosen point is then the next non-adjacent point, and so on. By default, the smallest \code{delta} 
is chosen to obtain a number of chosen points approximately equal to the specified \code{npts}.
Increasing \code{npts} or supplying a small \code{delta} will improve the accuracy of the fit (i.e. 
closer to the full lowess procedure) at the cost of running time. 

Robustification is performed using the magnitude of the residuals. Residuals greater than 6 times the 
median residual are assigned weights of zero. Otherwise, Tukey's biweight function is applied.
Weights are then used for weighted linear regression. Greater values of \code{iterations} will 
provide greater robustness.
}

\author{Aaron Lun}

\seealso{
\code{\link{lowess}}
}

\examples{
y <- rt(100,df=4)
x <- runif(100)
w <- runif(100)
out <- weightedLowess(x, y, w, span=0.7)
plot(x,y,cex=w)
o <- order(x)
lines(x[o],out$fitted[o],col="red")
}

\references{
Cleveland, W.S. (1979).
Robust Locally Weighted Regression and Smoothing Scatterplots.
\emph{Journal of the American Statistical Association} 74, 829-836.
}
