\name{03.ReadingData}
\alias{03.ReadingData}
\title{Topic: Reading Microarray Data from Files}

\description{
This help page gives an overview of LIMMA functions used to read data from files.
}

\section{Reading Target Information}{
The function \code{\link{readTargets}} is designed to help with organizing information about which RNA sample is hybridized to each channel on each array and which files store information for each array.
}

\section{Reading Intensity Data}{
The first step in a microarray data analysis is to read into R the intensity data for each array provided by an image analysis program.
This is done using the function \code{\link{read.maimages}}.

\code{\link{read.maimages}} optionally constructs quality weights for each spot using quality functions listed in \link{QualityWeights}.

If the data is two-color, then \code{read.maimages} produces an \code{RGList} object.
If the data is one-color (single channel) then an \code{EListRaw} object is produced.
In either case, \code{read.maimages} stores only the information required from each image analysis output file.
\code{\link{read.maimages}} uses utility functions \code{\link{removeExt}}, \code{\link{read.imagene}} and \code{\link{read.columns}}.
There are also a series of utility functions which read the header information from image output files including \code{\link{readGPRHeader}}, \code{\link{readImaGeneHeader}} and \code{\link{readGenericHeader}}.

\code{\link{read.ilmn}} reads probe or gene summary profile files from Illumina BeadChips,
and produces an \code{ElistRaw} object.

\code{\link{read.idat}} reads Illumina files in IDAT format, and produces an \code{EListRaw} object.
\code{\link{detectionPValues}} can be used to add detection p-values.

The function \link{as.MAList} can be used to convert a \code{marrayNorm} object to an \code{MAList} object if the data was read and normalized using the marray and marrayNorm packages.
}

\section{Reading the Gene List}{

Most image analysis software programs provide gene IDs as part of the intensity output files, for example GenePix, Imagene and the Stanford Microarray Database do this.
In other cases the probe ID and annotation information may be in a separate file.
The most common format for the probe annotation file is the GenePix Array List (GAL) file format.
The function \code{\link{readGAL}} reads information from a GAL file and produces a data frame with standard column names.

The function \code{\link{getLayout}} extracts from the GAL-file data frame the print layout information for a spotted array.
The functions \code{\link{gridr}}, \code{\link{gridc}}, \code{\link{spotr}} and \code{\link{spotc}} use the extracted layout to compute grid positions and spot positions within each grid for each spot.
The function \code{\link{printorder}} calculates the printorder, plate number and plate row and column position for each spot given information about the printing process.
The utility function \code{\link{getSpacing}} converts character strings specifying spacings of duplicate spots to numeric values.

The Australian Genome Research Facility in Australia often produces GAL files with composite probe IDs or names consisting of multiple strings separated by a delimiter.
These can be separated into name and annotation information using \code{\link{strsplit2}}.

If each probe is printed more than once of the arrays in a regular pattern, then \code{\link{uniquegenelist}} will remove duplicate names from the gal-file or gene list.
}

\section{Identifying Control Spots}{

The functions \code{\link{readSpotTypes}} and \code{\link{controlStatus}} assist with separating control spots from ordinary genes in the analysis and data exploration.
}

\section{Manipulating Data Objects}{

\code{\link[limma:cbind]{cbind}}, \code{\link[limma:cbind]{rbind}}, \code{\link[limma:merge]{merge}} allow different \code{RGList} or \code{MAList} objects to be combined.
\code{cbind} combines data from different arrays assuming the layout of the arrays to be the same.
\code{merge} can combine data even when the order of the probes on the arrays has changed.
\code{merge} uses utility function \code{\link{makeUnique}}.
}

\author{Gordon Smyth}

\seealso{
\link{01.Introduction},
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
\keyword{documentation}
