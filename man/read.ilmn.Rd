\name{read.ilmn}
\alias{read.ilmn}
\title{Read Illumina Expression Data}
\description{Read Illumina summary probe profile files and summary control probe profile files}
\usage{
read.ilmn(files=NULL, ctrlfiles=NULL, path=NULL, ctrlpath=NULL, probeid="Probe",
          annotation=c("TargetID", "SYMBOL"), expr="AVG_Signal",
          other.columns="Detection", sep="\t", quote="\"", verbose=TRUE, \dots)
}
\arguments{
  \item{files}{character vector giving the names of the summary probe profile files.}
  \item{ctrlfiles}{character vector giving the names of the summary control probe profile files.}
  \item{path}{character string giving the directory containing the summary probe profile files. Default is the current working directory.}
  \item{ctrlpath}{character string giving the directory containing the summary control probe profile files. Default is the same directory as for the probe profile files.}
  \item{probeid}{character string giving the name of the probe identifier column.}
  \item{annotation}{character vector giving possible column names for probe annotation.}
  \item{expr}{character string giving a keyword identifying the expression intensity columns. Any input column with column name containing this key will be read as containing intensity values.}
  \item{other.columns}{character vector giving keywords sufficient to identify any extra data columns that should be read in, such as "Detection", "Avg_NBEADS", "BEAD_STDEV" etc. The default of \code{Detection} is usually sufficient to identify the columns containing detection p-values.}
  \item{sep}{the field separator character.}
  \item{quote}{character string of characters to be treated as quote marks.}
  \item{verbose}{logical, \code{TRUE} to report names of profile files being read.}
  \item{\dots}{any other parameters are passed on to \code{\link{read.columns}}.}
  }
\details{
Illumina BeadStudio ouputs probe intensities (regular probe intensities) and control probe intensities to summary probe profile files (containing regular probes) and summary control probe profile files, respectively.
If both \code{files} and \code{ctrlfiles} are not \code{NULL}, this function will combine the data read from the two file types and save them to an \code{\link{EListRaw-class}} object.
If one of them is \code{NULL}, then only the required data are read in.

Probe types are indicated in the \code{Status} column of \code{genes}, a component of the returned \code{\link{EListRaw-class}} object.
There are totally seven types of control probes including \code{negative}, \code{biotin}, \code{labeling}, \code{cy3_hyb}, \code{housekeeping}, \code{high_stringency_hyb} or \code{low_stringency_hyb}.
Regular probes have the probe type \code{regular}.
The \code{Status} column will not be created if \code{ctrlfiles} is \code{NULL}.

To read in columns other than \code{probeid}, \code{annotation} and \code{expr}, users needs to specify keywords in \code{other.columns}.
One keyword corresponds to one type of columns.
Examples of keywords are "Detection", "Avg_NBEADS", "BEAD_STDEV" etc.
}

\value{ 
An \code{\link{EListRaw-class}} object with the following components:
\item{E}{numeric matrix of intensities.}
\item{genes}{data.frame of probe annotation. Contains any columns specified by \code{annotation} that are found in the input files.}
\item{other}{a list of matrices corresponding to any \code{other.columns} found in the input files.}
}

\author{Wei Shi and Gordon K Smyth}

\seealso{ 

\code{\link{read.ilmn.targets}} reads in Illumina expression data using the file information extracted from a target data frame which is often created by the \code{\link{readTargets}} function.

\code{\link{neqc}} performs normexp by control background correction, log transformation and quantile between-array normalization for Illumina expression data.

\code{\link{normexp.fit.control}} estimates the parameters of the normal+exponential convolution model with the help of negative control probes.

\code{\link{propexpr}} estimates the proportion of expressed probes in a microarray.
}

\examples{
\dontrun{
x <- read.ilmn(files="sample probe profile.txt",
               ctrlfiles="control probe profile.txt")
}
# See neqc and beadCountWeights for other examples using read.ilmn
}

\keyword{reading data}
\concept{illumina microarrays}
\concept{microarray data file}
