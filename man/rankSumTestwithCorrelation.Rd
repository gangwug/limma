\name{rankSumTestWithCorrelation}
\alias{rankSumTestWithCorrelation}

\title{Two Sample Wilcoxon-Mann-Whitney Rank Sum Test Allowing For Correlation}

\description{
A extension of the well-known rank-based test, but allowing for correlations between cases.
}

\usage{
rankSumTestWithCorrelation(index, statistics, correlation=0, df=Inf)
}

\arguments{
  \item{index}{any index vector such that \code{statistics[index]} contains the values of the statistic for the test group.}
  \item{statistics}{numeric vector giving values of the test statistic.}
  \item{correlation}{numeric scalar, average correlation between cases in the test group.  Cases in the second group are assumed independent of each other and other the first group.} 
  \item{df}{degrees of freedom which the correlation has been estimated.} 
}

\value{
Numeric vector of length 2 containing the \code{left.tail} and \code{right.tail} p-values.
}

\details{
This function implements a correlation-adjusted version of the Wilcoxon-Mann-Whitney test proposed by Wu and Smyth (2012).
It tests whether the mean rank of statistics in the test group is greater or less than the mean rank of the remaining statistic values.

When the correlation (or variance inflation factor) is zero, the function performs the usual two-sample Wilcoxon-Mann-Whitney rank sum test.
The Wilcoxon-Mann-Whitney test is implemented following the formulas given in Zar (1999) Section 8.10, including corrections for ties and for continuity.

The test allows for the possibility that cases in the test group may be more highly correlated on average than cases not in the group.
When the correlation is non-zero, the variance of the rank-sum statistic is computing using a formula derived from equation (4.5) of Barry et al (2008).
When the correlation is positive, the variance is increased and test will become more conservative.
}

\seealso{
\code{\link{wilcox.test}} performs the usual Wilcoxon-Mann-Whitney test assuming independence.

An overview of tests in limma is given in \link{08.Tests}.
}

\author{Gordon Smyth and Di Wu}

\references{
Barry, W.T., Nobel, A.B., and Wright, F.A. (2008).
A statistical framework for testing functional categories in microarray data.
\emph{Annals of Applied Statistics} 2, 286-315.

Wu, D, and Smyth, GK (2012). Camera: a competitive gene set test accounting for inter-gene correlation.
\emph{Nucleic Acids Research} 40, e133.
\url{http://nar.oxfordjournals.org/content/40/17/e133}

Zar, JH (1999). \emph{Biostatistical Analysis 4th Edition}.
Prentice-Hall International, Upper Saddle River, New Jersey.
}

\examples{
stat <- rnorm(100)
index <- 1:10
stat[index] <- stat[1:10]+1

rankSumTestWithCorrelation(index, stat)
rankSumTestWithCorrelation(index, stat, correlation=0.1)

group <- rep(1,100)
group[index] <- 2
group <- factor(group)
wilcox.test(stat ~ group)
}

\keyword{htest}

