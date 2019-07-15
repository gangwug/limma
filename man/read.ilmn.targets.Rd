\name{read.ilmn.targets}
\alias{read.ilmn.targets}
\title{Read Illumina Data from a Target Dataframe}
\description{Read Illumina data from a target dataframe}
\usage{
read.ilmn.targets(targets, \dots)
}
\arguments{
  \item{targets}{ data frame including names of profile files.}
  \item{\dots}{ any other parameters are passed on to \code{\link{read.ilmn}}.}
  }
\details{
\code{targets} is often created by calling the function \code{\link{readTargets}}.
Rows in \code{targets} are arrays and columns contain related array or RNA sample information.

At least one of the two columns called \code{files} and/or \code{ctrlfiles} should be present in \code{targets}, which includes names of summary probe profile files and names of summary control probe profile files respectively.
This function calls \code{\link{read.ilmn}} to read in the data.
}

\value{ 
An \code{\link{EListRaw-class}} object. See return value of the function \code{\link{read.ilmn}} for details.
}

\author{Wei Shi}

\seealso{ 
\code{\link{read.ilmn}}
}

\keyword{reading data}
\concept{illumina microarrays}
