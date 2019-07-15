\name{classifyTestsF}
\alias{classifyTestsF}
\title{Genewise Nested F-Tests}
\description{
For each gene, classify a series of related t-statistics as significantly up or down using nested F-tests.
}
\usage{
classifyTestsF(object, cor.matrix = NULL, df = Inf, p.value = 0.01, fstat.only = FALSE)
}
\arguments{
  \item{object}{numeric matrix of t-statistics or an \code{MArrayLM} object from which the t-statistics may be extracted.}
  \item{cor.matrix}{covariance matrix of each row of t-statistics. Will be extracted automatically from an \code{MArrayLM} object but otherwise defaults to the identity matrix.}
  \item{df}{numeric vector giving the degrees of freedom for the t-statistics.
  May have length 1 or length equal to the number of rows of \code{tstat}.
  Will be extracted automatically from an \code{MArrayLM} object but otherwise default to \code{Inf}.}
  \item{p.value}{numeric value between 0 and 1 giving the desired size of the test.}
  \item{fstat.only}{logical, if \code{TRUE} then return the overall F-statistic as for \code{FStat} instead of classifying the test results.}
}
\value{
If \code{fstat.only=FALSE}, then an object of class \code{\link[=TestResults-class]{TestResults}} is returned.
This is essentially a numeric matrix with elements \code{-1}, \code{0} or \code{1} depending on whether each t-statistic is classified as significantly negative, not significant or significantly positive respectively.

If \code{fstat.only=TRUE}, then a numeric vector of F-statistics is returned with attributes \code{df1} and \code{df2} giving the corresponding degrees of freedom.
}
\details{
\code{classifyTestsF} implements the \code{"nestedF"} multiple testing option offered by \code{\link{decideTests}}.
Users should generally use \code{decideTests} rather than calling \code{classifyTests} directly because, by itself, \code{classifyTests} does not incorporate any multiple testing adjustment across genes.
Instead it simply adjusts for multiple testing across contrasts for each gene individually.

\code{classifyTestsF} uses a nested F-test approach giving particular attention to correctly classifying genes that have two or more significant t-statistics, i.e., which are differentially expressed in two or more conditions.
For each row of \code{tstat}, the overall F-statistics is constructed from the t-statistics as for \code{FStat}.
At least one constrast will be classified as significant if and only if the overall F-statistic is significant.
If the overall F-statistic is significant, then the function makes a best choice as to which t-statistics contributed to this result.
The methodology is based on the principle that any t-statistic should be called significant if the F-test is still significant for that row when all the larger t-statistics are set to the same absolute size as the t-statistic in question.

Compared to conventional multiple testing methods, the nested F-test approach achieves better consistency between related contrasts.
(For example, if B is judged to be different from C, then at least one of B or C should be different to A.)
The approach was first used by Michaud et al (2008).
The nested F-test approach provides \emph{weak} control of the family-wise error rate, i.e., it correctly controls the type I error rate of calling any contrast as significant if all the null hypotheses are true.
In other words, it provides error rate control at the overall F-test level but does not provide strict error rate control at the individual contrast level. 

Usually \code{object} is a limma linear model fitted object, from which a matrix of t-statistics can be extracted, but it can also be a numeric matrix of t-statistics.
In either case, rows correspond to genes and columns to coefficients or contrasts.
If \code{object} is a matrix, then it may be necessary to supply values for \code{cor.matrix} and \code{df}.
The \code{cor.matrix} is the same as the correlation matrix of the coefficients from which the t-statistics were calculated and \code{df} is the degrees of freedom of the t-statistics.
All statistics for the same gene must have the same degrees of freedom.

If \code{fstat.only=TRUE}, the \code{classifyTestsF} just returns the vector of overall F-statistics for each gene.
}
\references{
Michaud, J, Simpson, KM, Escher, R, Buchet-Poyau, K, Beissbarth, T, Carmichael, C, Ritchie, ME, Schutz, F, Cannon, P, Liu, M, Shen, X, Ito, Y, Raskind, WH, Horwitz, MS, Osato, M, Turner, DR, Speed, TP, Kavallaris, M, Smyth, GK, and Scott, HS (2008).
Integrative analysis of RUNX1 downstream pathways and target genes.
\emph{BMC Genomics} 9, 363. 
}
\seealso{
An overview of multiple testing functions is given in \link{08.Tests}.
}
\author{Gordon Smyth}
\examples{
TStat <- matrix(c(0,10,0, 0,5,0, -4,-4,4, 2,2,2), 4, 3, byrow=TRUE)
colnames(TStat) <- paste0("Contrast",1:3)
rownames(TStat) <- paste0("Gene",1:4)
classifyTestsF(TStat, df=20)
FStat <- classifyTestsF(TStat, df=20, fstat.only=TRUE)
P <- pf(FStat, df1=attr(FStat,"df1"), df2=attr(FStat,"df2"), lower.tail=FALSE)
data.frame(F.Statistic=FStat,P.Value=P)
}
\keyword{htest}
