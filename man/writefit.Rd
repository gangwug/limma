\name{write.fit}
\alias{write.fit}

\title{Write MArrayLM Object to a File}

\description{
Write a microarray linear model fit to a file.
}

\usage{
write.fit(fit, results = NULL, file, digits = NULL,
          adjust = "none", method = "separate", F.adjust = "none",
          quote = FALSE, sep = "\t", row.names = TRUE, \dots)
}

\arguments{
\item{fit}{object of class \code{MArrayLM} containing the results of a linear model fit}
\item{results}{object of class \code{TestResults}}
\item{file}{character string giving name of file}
\item{digits}{integer indicating rounding precision for output values. If \code{NULL}, then no rounding is done.}
\item{adjust}{character string specifying multiple-testing adjustment method for the t-statistic P-values, e.g., \code{"BH"}. See \code{\link[stats]{p.adjust}} for the available options. If \code{NULL} or \code{"none"} then the P-values are not adjusted.}
\item{method}{character string, should the P-value adjustment be \code{"global"} or \code{"separate"} for each contrast.}
\item{F.adjust}{character string specifying adjustment method for the F-statistic P-values.}
\item{quote}{logical. If \code{TRUE}, any character or factor columns will be surrounded by double quotes.}
\item{sep}{the field separator string. Values in the output file will be separated by this string.}
\item{row.names}{logical, whether to include row names in output file.}
\item{\dots}{other arguments are passed to \code{write.table}}
}

\value{
No value is produced but a file is written to the current working directory.
}

\details{
This function writes a tab-delimited text file containing for each gene (1) the average log2-intensity, (2) the coefficients or contrasts (log2-fold-changes), (3) moderated t-statistics, (4) t-statistic P-values, (5) F-statistic if available, (6) F-statistic P-values if available, (7) classification if available and (8) gene names and annotation.
}

\author{Gordon Smyth}
\seealso{
\code{\link{write.table}} in the base library.

An overview of linear model functions in limma is given by \link{06.LinearModels}.
}
\keyword{IO}
