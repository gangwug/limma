\name{read.maimages}
\alias{read.maimages}
\alias{read.imagene}
\title{Read RGList or EListRaw from Image Analysis Output Files}
\description{
Reads an RGList from a set of two-color microarray image analysis output files,
or an EListRaw from a set of one-color files.
}
\usage{
read.maimages(files=NULL, source="generic", path=NULL, ext=NULL, names=NULL,
              columns=NULL, other.columns=NULL, annotation=NULL, green.only=FALSE,
              wt.fun=NULL, verbose=TRUE, sep="\t", quote=NULL, \dots)
read.imagene(files, path=NULL, ext=NULL, names=NULL, columns=NULL, other.columns=NULL,
             wt.fun=NULL, verbose=TRUE, sep="\t", quote="\"", \dots)
}
\arguments{
  \item{files}{character vector giving the names of the files containing image analysis output or, for Imagene data, a character matrix of names of files.
  Alternatively, it can be a data.frame containing a column called \code{FileName}.
  If omitted, then all files with extension \code{ext} in the specified directory will be read in alphabetical order.}
  \item{source}{character string specifying the image analysis program which produced the output files.  Choices are \code{"generic"}, \code{"agilent"}, \code{"agilent.median"}, \code{"agilent.mean"}, \code{"arrayvision"}, \code{"arrayvision.ARM"}, \code{"arrayvision.MTM"}, \code{"bluefuse"}, \code{"genepix"}, \code{"genepix.custom"}, \code{"genepix.median"}, \code{"imagene"}, \code{"imagene9"}, \code{"quantarray"}, \code{"scanarrayexpress"}, \code{"smd.old"}, \code{"smd"}, \code{"spot"} or \code{"spot.close.open"}.}
  \item{path}{character string giving the directory containing the files.
  The default is the current working directory.}
  \item{ext}{character string giving optional extension to be added to each file name}
  \item{names}{character vector of unique names to be associated with each array as column name.
  Can be supplied as \code{files$Label} if \code{files} is a data.frame.
  Defaults to \code{removeExt(files)}.}
  \item{columns}{list, or named character vector.
  For two color data, this should have fields \code{R}, \code{G}, \code{Rb} and \code{Gb} giving the column names to be used for red and green foreground and background or, in the case of Imagene data, a list with fields \code{f} and \code{b}.
  For single channel data, the fields are usually \code{E} and \code{Eb}.
  This argument is optional if \code{source} is specified, otherwise it is required.}
  \item{other.columns}{character vector of names of other columns to be read containing spot-specific information}
  \item{annotation}{character vector of names of columns containing annotation information about the probes}
  \item{green.only}{logical, for use with \code{source}, should the green (Cy3) channel only be read, or are both red and green required?}
  \item{wt.fun}{function to calculate spot quality weights}
  \item{verbose}{logical, \code{TRUE} to report each time a file is read}
  \item{sep}{the field separator character}
  \item{quote}{character string of characters to be treated as quote marks}
  \item{\dots}{any other arguments are passed to \code{read.table}}
}
\details{
These are the main data input functions for the LIMMA package.
\code{read.maimages} reads either single channel or two-color microarray intensity data from text files.
\code{read.imagene} is specifically for two-color ImaGene intensity data created by ImaGene versions 1 through 8, and is called by \code{read.maimages} to read such data.

\code{read.maimages} is designed to read data from any microarray platform except for Illumina BeadChips, which are read by \code{\link{read.ilmn}}, and Affymetrix GeneChip data, which is best read and pre-processed by specialist packages designed for that platform.

\code{read.maimages} extracts the foreground and background intensities from a series of files, produced by an image analysis program, and assembles them into the components of one list.
The image analysis programs Agilent Feature Extraction, ArrayVision, BlueFuse, GenePix, ImaGene, QuantArray (Version 3 or later), Stanford Microarray Database (SMD) and SPOT are supported explicitly.
Almost all these programs write the intensity data for each microarray to one file.
The exception is ImaGene, early versions of which wrote the red and green channels of each microarray to different files.
Data from some other image analysis programs not mentioned above can be read if the appropriate column names containing the foreground and background intensities are specified using the \code{columns} argument.
(Reading custom columns will work provided the column names are unique and there are no rows in the file after the last line of data.
Header lines are ok.)

For Agilent files, two possible foreground estimators are supported: \code{source="agilent.median"} use median foreground while \code{source="agilent.mean"} uses mean foreground.
Background estimates are always medians.
The use of \code{source="agilent"} defaults to \code{"agilent.median"}.
Note that this behavior is new from 9 March 2012.
Previously, in limma 3.11.16 or earlier, \code{"agilent"} had the same meaning as \code{"agilent.mean"}.

For GenePix files, two possible foreground estimators are supported as well as custom background: \code{source="genepix.median"} uses the median foreground estimates while \code{source="genepix.mean"} uses mean foreground estimates.
The use of \code{source="genepix"} defaults to \code{"genepix.mean"}.
Background estimates are always medians unless \code{source="genepix.custom"} is specified.
GenePix 6.0 and later supply some custom background options, notably morphological background.
If the GPR files have been written using a custom background, then \code{source="genepix.custom"} will cause it to be read and used. 

For SPOT files, two possible background estimators are supported:
\code{source="spot"} uses background intensities estimated from the morphological opening algorithm.
If \code{source="spot.close.open"} then background intensities are estimated from morphological closing followed by opening.

ArrayVision reports spot intensities in a number of different ways.
\code{read.maimages} caters for ArrayVision's Artifact-removed (ARM) density values using \code{source="arrayvision.ARM"} or for
Median-based Trimmed Mean (MTM) density values with \code{"arrayvision.MTM"}.
ArrayVision users may find it useful to read the top two lines of their data file to check which version of density values they have.

SMD data should consist of raw data files from the database, in tab-delimited text form.
There are two possible sets of column names depending on whether the data was entered into the database before or after September 2003.
\code{source="smd.old"} indicates that column headings in use prior to September 2003 should be used.

Intensity data from ImaGene versions 1 to 8 (\code{source="imagene"}) is different from other image analysis programs in that the read and green channels were written to separate files.
\code{read.maimages} handles the special behaviour of the early ImaGene versions by requiring that the argument \code{files} should be a matrix with two columns instead of a vector.
The first column should contain the names of the files containing green channel (cy3) data and the second column should contain names of files containing red channel (cy5) data.
Alternately, \code{files} can be entered as a vector of even length instead of a matrix.
In that case, each consecutive pair of file names is assumed to contain the green (cy3) and red (cy5) intensities respectively from the same array.
The function \code{read.imagene} is called by \code{read.maimages} when \code{source="imagene"}, so \code{read.imagene} does not need to be called directly by users.

ImaGene version~9 (\code{source="imagene9"}) reverts to the same behavior as the other image analysis programs.
For ImaGene~9, \code{files} is a vector of length equal to the number of microarrays, same as for other image analysis programs.

Spot quality weights may be extracted from the image analysis files using a weight function wt.fun.
\code{wt.fun} may be any user-supplied function which accepts a data.frame argument and returns a vector of non-negative weights.
The columns of the data.frame are as in the image analysis output files.
There is one restriction, which is that the column names should be refered to in full form in the weight function, i.e., do not rely on name expansion for partial matches when refering to the names of the columns.
See \code{\link{QualityWeights}} for suggested weight functions.

The argument \code{other.columns} allows arbitrary columns of the image analysis output files to be preserved in the data object.
These become matrices in the component \code{other} component.
For ImaGene data, the other column headings should be prefixed with \code{"R "} or \code{"G "} as appropriate.
}

\section{Warnings}{
All image analysis files being read are assumed to contain data for the same genelist in the same order.
No checking is done to confirm that this is true.
Probe annotation information is read from the first file only.
}

\value{
For one-color data, an \code{\link[limma:EList]{EListRaw}} object.
For two-color data, an \code{\link[limma:rglist]{RGList}} object containing the components
  \item{R}{matrix containing the red channel foreground intensities for each spot for each array.}
  \item{Rb}{matrix containing the red channel background intensities for each spot for each array.}
  \item{G}{matrix containing the green channel foreground intensities for each spot for each array.}
  \item{Gb}{matrix containing the green channel background intensities for each spot for each array.}
  \item{weights}{spot quality weights, if \code{wt.fun} is given}
  \item{other}{list containing matrices corresponding to \code{other.columns} if given}
  \item{genes}{data frame containing annotation information about the probes, for example gene names and IDs and spatial positions on the array, currently set only if \code{source} is \code{"agilent"}, \code{"genepix"} or \code{source="imagene"} or if the \code{annotation} argument is set}
  \item{targets}{data frame with column \code{FileName} giving the names of the files read.  If \code{files} was a data.frame on input, then the whole data.frame is stored here on output.}
  \item{source}{character string giving the image analysis program name}
  \item{printer}{list of class \code{\link[=PrintLayout-class]{PrintLayout}}, currently set only if \code{source="imagene"}}
}

\author{Gordon Smyth, with speed improvements suggested by Marcus Davy}

\references{
Ritchie, ME, Phipson, B, Wu, D, Hu, Y, Law, CW, Shi, W, and Smyth, GK (2015).
limma powers differential expression analyses for RNA-sequencing and microarray studies.
\emph{Nucleic Acids Research} 43, e47.
\url{http://nar.oxfordjournals.org/content/43/7/e47}

Web pages for the image analysis software packages mentioned here are listed at \url{http://www.statsci.org/micrarra/image.html}
}
\seealso{
\code{read.maimages} uses \code{\link{read.columns}} for efficient reading of text files.
As far as possible, it is has similar behavior to \code{\link[utils]{read.table}} in the base package.

\code{\link{read.ilmn}} reads probe or gene summary profile files from Illumina BeadChips.

An overview of LIMMA functions for reading data is given in \link{03.ReadingData}.
}
\examples{
#  Read all .gpr files from current working directory
#  and give weight 0.1 to spots with negative flags

\dontrun{files <- dir(pattern="*\\\\.gpr$")
RG <- read.maimages(files,"genepix",wt.fun=wtflags(0.1))}

#  Read all .spot files from current working director and down-weight
#  spots smaller or larger than 150 pixels

\dontrun{files <- dir(pattern="*\\\\.spot$")
RG <- read.maimages(files,"spot",wt.fun=wtarea(150))}
}

\keyword{reading data}
\concept{microarray data file}
