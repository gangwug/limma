#	CLASSES.R

setClass("RGList",
#  Class to hold initial read-in two-color data
representation("list")
)

setClass("MAList",
#  Class to hold normalized, rotated two-color data
representation("list")
)

setClass("EListRaw",
#  Class to hold one channel data on raw scale
representation("list")
)

setClass("EList",
#  Class to hold one channel data on log scale
representation("list")
)

setClass("MArrayLM",
#  Linear model fit
representation("list")
)

printHead <- function(x)
#  Print leading 5 elements or rows of atomic object
#  Gordon Smyth
#  May 2003.  Last modified 4 March 2017.
{
	if(is.atomic(x)) {
		d <- dim(x)
		if(length(d)<2) which <- "OneD"
		if(length(d)==2) which <- "TwoD"
		if(length(d)>2) which <- "Array"
	} else {
		if(inherits(x,"data.frame")) {
			d <- dim(x)
			which <- "TwoD"
		} else {
			if(is.call(x))
				which <- "Call"
			else {
				if(is.recursive(x))
					which <- "Recursive"
				else
					which <- "Other"
			}
		}
	}
	switch(which,
	OneD={
		n <- length(x)
		if(n > 20) {
			print(x[1:5])
			cat(n-5,"more elements ...\n")
		} else
			print(x)
	},
	TwoD={
		n <- d[1]
		if(n > 10) {
			print(x[1:5,,drop=FALSE])
			cat(n-5,"more rows ...\n")
		} else
			print(x)
	},
	Array={
		n <- d[1]
		if(n > 10) {
			dn <- dimnames(x)
			dim(x) <- c(d[1],prod(d[-1]))
			x <- x[1:5,,drop=FALSE]
			dim(x) <- c(5,d[-1])
			if(!is.null(dn[[1]])) dn[[1]] <- dn[[1]][1:5]
			dimnames(x) <- dn
			print(x)
			cat(n-5,"more rows ...\n")
		} else
			print(x)
	},
	Recursive={
		n <- length(x)
		if(n) {
			i <- names(x)
			if(is.null(i)) i <- seq_len(n)
			for (what in i) {
				y <- x[[what]]
				cat("$",what,"\n",sep="")
				Recall(y)
				cat("\n")
			}
		}
	},
	Call=,Other=print(x)
	)
}

setClass("LargeDataObject")
setIs("RGList","LargeDataObject")
setIs("MAList","LargeDataObject")
setIs("EListRaw","LargeDataObject")
setIs("EList","LargeDataObject")
setIs("MArrayLM","LargeDataObject")

setMethod("show","LargeDataObject",
#  Print and show method large data objects
#  Gordon Smyth
#  May 2003
function(object) {
	cat("An object of class \"",class(object),"\"\n",sep="")
	for (what in names(object)) {
		x <- object[[what]]
		cat("$",what,"\n",sep="")
		printHead(x)
		cat("\n")
	}
	for (what in setdiff(slotNames(object),".Data")) {
		x <- slot(object,what)
		if(length(x) > 0) {
			cat("@",what,"\n",sep="")
			printHead(x)
			cat("\n")
		}
	}
})

dim.RGList <- function(x) if(is.null(x$R)) c(0,0) else dim(as.matrix(x$R))
dim.MAList <- function(x) if(is.null(x$M)) c(0,0) else dim(as.matrix(x$M))
dim.EListRaw <- dim.EList <- function(x) if(is.null(x$E)) c(0,0) else dim(as.matrix(x$E))
dim.MArrayLM <- function(x) if(is.null(x$coefficients)) c(0,0) else dim(as.matrix(x$coefficients))
#length.RGList <- length.MAList <- length.EListRaw <- length.EList <- length.MArrayLM <- function(x) prod(dim(x))

dimnames.RGList <- function(x) dimnames(x$R)
dimnames.MAList <- function(x) dimnames(x$M)
dimnames.EListRaw <- dimnames.EList <- function(x) dimnames(x$E)
dimnames.MArrayLM <- function(x) dimnames(x$coefficients)
.setdimnames <- function(x, value)
#  Dimension names for RGList-like objects
#  Gordon Smyth
#  17 Dec 2005. Last modified 23 March 2009.
{
	exprmatrices <- c("R","G","Rb","Gb","M","A","E","Eb","weights")
	for (a in exprmatrices) if(!is.null(x[[a]])) dimnames(x[[a]]) <- value
	for(a in names(x$other)) dimnames(x$other[[a]]) <- value
	if(!is.null(x$targets)) row.names(x$targets) <- value[[2]]
	if(!is.null(x$design)) rownames(x$design) <- value[[2]]
	x
}
#assign("dimnames<-.RGList",.setdimnames)
#assign("dimnames<-.MAList",.setdimnames)
"dimnames<-.RGList" <- .setdimnames
"dimnames<-.MAList" <- .setdimnames
"dimnames<-.EListRaw" <- "dimnames<-.EList" <- .setdimnames

summary.MArrayLM <- summary.MAList <- summary.RGList <- summary.EListRaw <- summary.EList <- function(object,...) summary(unclass(object))

as.MAList <- function(object) {
#	Convert marrayNorm object to MAList
#	Gordon Smyth
#	20 Sep 2003.  Last modified 20 Dec 2003.

	MA <- new("MAList")
	ifposlen <- function(x) if(length(x)) return(x) else return(NULL)
	MA$A <- ifposlen(object@maA)
	MA$M <- ifposlen(object@maM)
	MA$weights <- ifposlen(object@maW)
	MA$printer$ngrid.r <- ifposlen(object@maLayout@maNgr)
	MA$printer$ngrid.c <- ifposlen(object@maLayout@maNgc)
	MA$printer$nspot.r <- ifposlen(object@maLayout@maNsr)
	MA$printer$nspot.c <- ifposlen(object@maLayout@maNsc)
	MA$printer$notes <- ifposlen(object@maLayout@maNotes)
	MA$genes <- ifposlen(object@maGnames@maInfo)
	MA$genes$Labels <- ifposlen(object@maGnames@maLabels)
	attr(MA$genes,"notes") <- ifposlen(object@maGnames@maNotes)
	MA$genes$Sub <- ifposlen(object@maLayout@maSub)
	MA$genes$Plate <- ifposlen(object@maLayout@maPlate)
	MA$genes$Controls <- ifposlen(object@maLayout@maControls)
	MA$targets <- ifposlen(object@maTargets@maInfo)
	MA$targets$Labels <- ifposlen(object@maTargets@maLabels)
	MA$notes <- ifposlen(object@maNotes)
	MA$maNormCall <- ifposlen(object@maNormCall)
	MA
} 

#  Gordon Smyth, 28 Oct 2004, 23 March 2009
as.matrix.RGList <- function(x,...) normalizeWithinArrays(x,method="none")$M
as.matrix.MAList <- function(x,...) as.matrix(x$M)
as.matrix.EListRaw <- as.matrix.EList <- function(x,...) as.matrix(x$E)
as.matrix.MArrayLM <- function(x,...) x$coefficients
as.matrix.marrayNorm <- function(x,...) x@maM
#  13 July 2006
as.matrix.PLMset <- function(x,...) x@chip.coefs
#  19 Dec 2006, 18 May 2007
as.matrix.ExpressionSet <- as.matrix.LumiBatch <- function(x,...) env=x@assayData[["exprs"]]
#  16 Sep 2007
as.matrix.vsn <- function(x,...) x@hx

as.data.frame.MAList <- function(x, row.names = NULL, optional = FALSE, ...)
#	Convert MAList object to data.frame
#	Gordon Smyth
#	10 Dec 2009.
{
	if(is.null(row.names) && !is.null(rownames(x$M))) row.names <- makeUnique(rownames(x$M))
	if(is.null(x$genes)) {
		y <- as.data.frame(x$M,row.names=row.names,check.names=FALSE)
	} else {
		if(is.vector(x$genes)) x$genes <- data.frame(ID=x$genes,stringsAsFactors=FALSE)
		y <- data.frame(x$genes,x$M,row.names=row.names,check.names=FALSE,stringsAsFactors=FALSE)
	}
	y
}

as.data.frame.EList <- as.data.frame.EListRaw <- function(x, row.names = NULL, optional = FALSE, ...)
#	Convert EList object to data.frame
#	Gordon Smyth
#	11 Dec 2009.
{
	if(is.null(row.names) && !is.null(rownames(x$E))) row.names <- makeUnique(rownames(x$E))
	if(is.null(x$genes)) {
		y <- as.data.frame(x$E,row.names=row.names,check.names=FALSE)
	} else {
		if(is.vector(x$genes)) x$genes <- data.frame(ID=x$genes,stringsAsFactors=FALSE)
		y <- data.frame(x$genes,x$E,row.names=row.names,check.names=FALSE,stringsAsFactors=FALSE)
	}
	y
}

as.data.frame.MArrayLM <- function(x, row.names = NULL, optional = FALSE, ...)
#	Convert MArrayLM object to data.frame
#	Gordon Smyth
#	6 April 2005.  Last modified 13 Jan 2006.
{
	x <- unclass(x)
	if(is.null(x$coefficients)) {
		warning("NULL coefficients, returning empty data.frame")
		return(data.frame())
	}
	cn <- names(x)
	nprobes <- NROW(x$coefficients)
	ncoef <- NCOL(x$coefficients)
	include.comp <- cn[unlist(lapply(x,NROW))==nprobes]
	other.comp <- setdiff(names(x),include.comp)
	if(length(other.comp)) for (a in other.comp) x[[a]] <- NULL
#	coef.comp <- c("coefficients","stdev.unscaled","t","p.value","lods")
#	for (a in coef.comp) if(!is.null(x[[a]]) && NCOL(x[[a]])==1) colnames(x[[a]]) <- paste(a,colnames(x[[a]]),sep=".")
	if(ncoef==1) x <- lapply(x,drop)
	as.data.frame(x,row.names=row.names,optional=optional)
}
