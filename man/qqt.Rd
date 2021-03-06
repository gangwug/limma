\name{qqt}
\alias{qqt}
\alias{qqf}
\title{Student's t or Fisher's F Quantile-Quantile Plot}
\description{Plots the quantiles of a data sample against the theoretical quantiles of a Student's t distribution.}
\usage{qqt(y, df = Inf, ylim = range(y), main = "Student's t Q-Q Plot", 
    xlab = "Theoretical Quantiles", ylab = "Sample Quantiles", plot.it = TRUE, \dots) 
qqf(y, df1, df2, ylim=range(y), main= "F Distribution Q-Q Plot",
    xlab = "Theoretical Quantiles", ylab = "Sample Quantiles", plot.it = TRUE, \dots)
}
\arguments{
\item{y}{a numeric vector or array containing the data sample}
\item{df}{degrees of freedom for the t-distribution.  The default \code{df=Inf} represents the normal distribution.}
\item{df1}{numerator degrees of freedom for the F-distribution.}
\item{df2}{denominator degrees of freedom for the F-distribution.}
\item{ylim}{plotting range for \code{y}}
\item{main}{main title for the plot}
\item{xlab}{x-axis title for the plot}
\item{ylab}{y-axis title for the plot}
\item{plot.it}{whether or not to produce a plot}
\item{\dots}{other arguments to be passed to \code{plot}}
}
\value{A list is invisibly returned containing the values plotted in the QQ-plot:
\item{x}{theoretical quantiles of the t-distribution or F-distribution}
\item{y}{the data sample, same as input \code{y}}
}
\details{
This function is analogous to \code{qqnorm} for normal probability plots.
In fact \code{qqt(y,df=Inf)} is identical to \code{qqnorm(y)} in all respects except the default title on the plot.
}
\author{Gordon Smyth}
\seealso{\code{\link[stats]{qqnorm}}}
\examples{
#  See also the lmFit examples

y <- rt(50,df=4)
qqt(y,df=4)
abline(0,1)
}

\keyword{distributions}
