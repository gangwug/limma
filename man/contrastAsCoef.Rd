\name{contrastAsCoef}
\alias{contrastAsCoef}
\title{Reform a Design Matrix to that Contrasts Become Coefficients}
\description{
Reform a design matrix so that one or more coefficients from the new matrix correspond to specified contrasts of coefficients from the old matrix.
}
\usage{
contrastAsCoef(design, contrast=NULL, first=TRUE)
}
\arguments{
  \item{design}{numeric design matrix.}
  \item{contrast}{numeric matrix with rows corresponding to columns of the design matrix (coefficients) and columns containing contrasts. May be a vector if there is only one contrast.}
  \item{first}{logical, should coefficients corresponding to contrasts be the first columns (\code{TRUE}) or last columns (\code{FALSE}) of the output design matrix.}
}

\details{
If the contrasts contained in the columns of \code{contrast} are not linearly dependent, then superfluous columns are dropped until the remaining matrix has full column rank.
The number of retained contrasts is stored in \code{qr$rank} and the retained columns are given by \code{qr$pivot}.
}

\value{
A list with components
\item{design}{reformed design matrix}
\item{coef}{columns of design matrix which hold the meaningful coefficients}
\item{qr}{QR-decomposition of contrast matrix}
}

\seealso{
\code{\link[stats]{model.matrix}} in the stats package.

An overview of linear model functions in limma is given by \link{06.LinearModels}.
}

\author{Gordon Smyth}

\examples{
design <- cbind(1,c(0,0,1,1,0,0),c(0,0,0,0,1,1))
cont <- c(0,-1,1)
design2 <- contrastAsCoef(design, cont)$design

#  Original coef[3]-coef[2] becomes coef[1]
y <- rnorm(6)
fit1 <- lm(y~0+design)
fit2 <- lm(y~0+design2)
coef(fit1)
coef(fit1) %*% cont
coef(fit2)
}

\keyword{regression}
