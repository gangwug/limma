\name{read.idat}
\alias{read.idat}

\title{Read Illumina expression data directly from IDAT files}

\description{Read Illumina BeadArray data from IDAT and manifest (.bgx) files for gene expression platforms.}

\usage{
read.idat(idatfiles, bgxfile, dateinfo = FALSE, annotation = "Symbol",
          tolerance = 0L, verbose = TRUE)
}

\arguments{
  \item{idatfiles}{character vector specifying idat files to be read in.}
  \item{bgxfile}{character string specifying bead manifest file (.bgx) to be read in.}
  \item{dateinfo}{logical. Should date and software version information be read in?}
  \item{annotation}{character vector of annotation columns to be read from the manifest file.}
  \item{tolerance}{integer. The number of probe ID discrepancies allowed between the manifest and any of the IDAT files.} 
  \item{verbose}{logical. Should progress messages are sent to standard output?}
}

\details{
     Illumina's BeadScan/iScan software outputs probe intensities in IDAT
     format (encrypted XML files) and uses probe information stored in a platform specific manifest file (.bgx).
     These files can be processed using the low-level functions \code{readIDAT} and \code{readBGX} 
     from the \code{illuminaio} package (Smith et al. 2013).

     The \code{read.idat} function provides a convenient way to read these files
     into R and to store them in an \code{EListRaw-class} object.
     The function serves a similar purpose to \code{\link{read.ilmn}}, 
     which reads text files exported by Illumina's GenomeStudio software,
     but it reads the IDAT files directly without any need to convert them first to text.
 
     The function reads information on control probes as well for regular probes.
     Probe types are indicated in the \code{Status} column of the \code{genes} 
     component of the \code{EListRaw} object.

     The \code{annotation} argument specifies probe annotation columns to be extracted from the manifest file.
     The manifest typically contains the following columns:
 \code{"Species"}, \code{"Source"}, \code{"Search_Key"}, \code{"Transcript"},
 \code{"ILMN_Gene"}, \code{"Source_Reference_ID"}, \code{"RefSeq_ID"},
 \code{"Unigene_ID"}, \code{"Entrez_Gene_ID"}, \code{"GI"},
 \code{"Accession"}, \code{"Symbol"}, \code{"Protein_Product"},
 \code{"Probe_Id"}, \code{"Array_Address_Id"}, \code{"Probe_Type"},
 \code{"Probe_Start"}, \code{"Probe_Sequence"}, \code{"Chromosome"},
 \code{"Probe_Chr_Orientation"}, \code{"Probe_Coordinates"}, \code{"Cytoband"},
 \code{"Definition"}, \code{"Ontology_Component"}, \code{"Ontology_Process"},
 \code{"Ontology_Function"}, \code{"Synonyms"}, \code{"Obsolete_Probe_Id"}.
     Note that \code{"Probe_Id"} and \code{"Array_Address_Id"} are always extracted and
     do not need to included in the \code{annotation} argument.

     If more than \code{tolerance} probes in the manifest cannot be found in an IDAT file then the function will return an error.
}

\value{
  An \code{EListRaw} object with the following components:
  \item{E}{ numeric matrix of raw intensities.}
  \item{other$NumBeads}{ numeric matrix of same dimensions as \code{E} giving number of beads used for each intensity value.}
  \item{other$STDEV}{ numeric matrix of same dimensions as \code{E} giving bead-level standard deviation or standard error for each intensity value.}
  \item{genes}{ data.frame of probe annotation.
  This includes the \code{Probe_Id} and \code{Array_Address_Id} columns extracted from the manifest file,
  plus a \code{Status} column identifying control probes,
  plus any other columns specified by \code{annotation}.}
  \item{targets}{ data.frame of sample information.
  This includes the IDAT file names plus other columns if \code{dateinfo=TRUE}.}
}

\references{
Smith ML, Baggerly KA, Bengtsson H, Ritchie ME, Hansen KD (2013). 
illuminaio: An open source IDAT parsing tool. \emph{F1000 Research} 2, 264.
\url{http://f1000research.com/articles/2-264/}
}

\author{Matt Ritchie}

\seealso{
     \code{\link{read.ilmn}} imports gene expression data output by GenomeStudio.

     \code{\link{neqc}} performs normexp by control background correction, log
     transformation and quantile between-array normalization for
     Illumina expression data.

     \code{\link{propexpr}} estimates the proportion of expressed probes in a microarray.
     
     \code{\link{detectionPValues}} computes detection p-values from the negative controls.
}

\examples{
\dontrun{
idatfiles <- dir(pattern="idat")
bgxfile <- dir(pattern="bgx")
x <- read.idat(idatfiles, bgxfile)
x$other$Detection <- detectionPValues(x)
propexpr(data)
y <- neqc(data)
}
}

\keyword{reading data}
\concept{illumina microarrays}
\concept{microarray data file}
