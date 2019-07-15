\name{01.Introduction}
\docType{package}
\alias{01.Introduction}
\alias{limma}
\alias{limma-package}
\title{Introduction to the LIMMA Package}
\description{
LIMMA is a library for the analysis of gene expression microarray data, especially the use of linear models for analysing designed experiments and the assessment of differential expression.
LIMMA provides the ability to analyse comparisons between many RNA targets simultaneously in arbitrary complicated designed experiments.
Empirical Bayesian methods are used to provide stable results even when the number of arrays is small.
The linear model and differential expression functions apply to all gene expression technologies, including microarrays, RNA-seq and quantitative PCR.
}

\details{
There are three types of documentation available:
\enumerate{
\item
The \emph{LIMMA User's Guide} can be reached through the "User
Guides and Package Vignettes" links at the top of the LIMMA
contents page.  The function \code{\link{limmaUsersGuide}} gives
the file location of the User's Guide.\cr
\item
An overview of limma functions grouped by purpose is contained
in the numbered chapters at the foot of the LIMMA package index page,
of which this page is the first.\cr
\item
The LIMMA contents page gives an
alphabetical index of detailed help topics.\cr
}

The function \code{\link{changeLog}} displays the record of changes to the package.
}

\author{Gordon Smyth, with contributions from many colleagues}

\references{
Phipson, B, Lee, S, Majewski, IJ, Alexander, WS, and Smyth, GK (2016).
Robust hyperparameter estimation protects against hypervariable genes and improves power to detect differential expression.
\emph{Annals of Applied Statistics} 10, 946-963.
\url{http://projecteuclid.org/euclid.aoas/1469199900}

Ritchie, ME, Phipson, B, Wu, D, Hu, Y, Law, CW, Shi, W, and Smyth, GK (2015).
limma powers differential expression analyses for RNA-sequencing and microarray studies.
\emph{Nucleic Acids Research} 43, e47.
\url{http://nar.oxfordjournals.org/content/43/7/e47}

Law, CW, Chen, Y, Shi, W, and Smyth, GK (2014).
Voom: precision weights unlock linear model analysis tools for RNA-seq read counts.
\emph{Genome Biology} 15, R29.
\url{http://genomebiology.com/2014/15/2/R29}

Smyth, G. K. (2004). Linear models and empirical Bayes methods for assessing differential expression in microarray experiments.
\emph{Statistical Applications in Genetics and Molecular Biology}, Volume 3, Article 3.
\url{http://www.statsci.org/smyth/pubs/ebayes.pdf}
}

\seealso{
\link{02.Classes},
\link{03.ReadingData},
\link{04.Background},
\link{05.Normalization},
\link{06.LinearModels},
\link{07.SingleChannel},
\link{08.Tests},
\link{09.Diagnostics},
\link{10.GeneSetTests},
\link{11.RNAseq}
}

\keyword{package}
