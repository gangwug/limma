##  plotMD.R

plotMD <- function(object,...) UseMethod("plotMD")

plotMD.RGList <- function(object, column=1, array=NULL, xlab="A", ylab="M", main=colnames(object)[column], status=object$genes$Status, zero.weights=FALSE, ...)
#	Mean-difference plot with color coding for controls
#	Created by Gordon Smyth 7 April 2003 and James Wettenhall 27 June 2003.
#	Last modified 7 June 2015.
{
	if(!is.null(array)) column <- array
	object <- MA.RG(object[,column])
	plotMD.MAList(object=object,column=1,xlab=xlab,ylab=ylab,main=main,status=status,zero.weights=zero.weights,...)
}

plotMD.MAList <- function(object, column=1, array=NULL, xlab="A", ylab="M", main=colnames(object)[column], status=object$genes$Status, zero.weights=FALSE, ...)
#	Mean-difference plot with color coding for controls
#	Gordon Smyth 7 April 2003, James Wettenhall 27 June 2003.
#	Last modified 7 June 2015.
{
	if(!is.null(array)) column <- array
	A <- as.matrix(object$A)[,column]
	M <- as.matrix(object$M)[,column]
	if(!zero.weights && !is.null(object$weights)) {
		w <- as.matrix(object$weights)[,column]
		M[ is.na(w) | (w <= 0) ] <- NA
	}
	plotWithHighlights(x=A,y=M,xlab=xlab,ylab=ylab,main=main,status=status,...)
}

plotMD.MArrayLM <- function(object, column=ncol(object), coef=NULL, xlab="Average log-expression", ylab="log-fold-change", main=colnames(object)[column], status=object$genes$Status, zero.weights=FALSE, ...)
#	Mean-difference plot with color coding for controls
#	Gordon Smyth 7 April 2003, James Wettenhall 27 June 2003.
#	Last modified 24 June 2015.
{
	if(!is.null(coef)) column <- coef
	if(is.null(object$Amean)) stop("Amean component is absent.")
	logFC <- as.matrix(object$coefficients)[,column]
	if(!zero.weights && !is.null(object$weights)) {
		w <- as.matrix(object$weights)[,column]
		logFC[ is.na(w) | (w <= 0) ] <- NA
	}
	plotWithHighlights(x=object$Amean,y=logFC,xlab=xlab,ylab=ylab,main=main,status=status,...)
}

plotMD.EListRaw <- function(object, column=1, array=NULL, xlab="Average log-expression", ylab="Expression log-ratio (this sample vs others)", main=colnames(object)[column], status=object$genes$Status, zero.weights=FALSE, ...)
#	Mean-difference plot with color coding for controls
#	Gordon Smyth
#	Created 22 Oct 2014. Last modified 7 June 2015.
{
	if(!is.null(array)) column <- array
	if(!is.null(object$Eb)) object$E <- object$E-object$Eb
	object$E <- log2(object$E)
	object <- new("EList",unclass(object))
	plotMD(object, column=column, xlab=xlab, ylab=ylab, main=main, status=status, zero.weights=zero.weights, ...)
}

plotMD.EList <- function(object, column=1, array=NULL, xlab="Average log-expression", ylab="Expression log-ratio (this sample vs others)", main=colnames(object)[column], status=object$genes$Status, zero.weights=FALSE, ...)
#	Mean-difference plot with color coding for controls
#	Gordon Smyth 7 April 2003, James Wettenhall 27 June 2003.
#	Last modified 14 April 2015.
{
	if(!is.null(array)) column <- array
	E <- as.matrix(object$E)
	if(ncol(E) < 2L) stop("Need at least two columns")

#	Convert column to integer if not already
	j <- 1:ncol(E)
	names(j) <- colnames(E)
	column <- j[column[1]]

	AveOfOthers <- rowMeans(E[,-column,drop=FALSE],na.rm=TRUE)
	Diff <- E[,column]-AveOfOthers
	Mean <- (E[,column]+AveOfOthers)/2

	if(!zero.weights && !is.null(object$weights)) {
		w <- as.matrix(object$weights)[,column]
		Diff[ is.na(w) | (w <= 0) ] <- NA
	}

	plotWithHighlights(x=Mean,y=Diff,xlab=xlab,ylab=ylab,main=main,status=status,...)
}

plotMD.default <- function(object, column=1, xlab="Average log-expression", ylab="Expression log-ratio (this sample vs others)", main=colnames(object)[column], status=NULL, ...)
#	Mean-difference plot with color coding for controls
#	Gordon Smyth 7 April 2003, James Wettenhall 27 June 2003.
#	Last modified 14 April 2015.
{
#	Data is assumed to be single-channel
	object <- as.matrix(object)
	ncolumns <- ncol(object)
	if(ncolumns<2) stop("Need at least two columns")
	column <- as.integer(column[1L])
	Ave <- rowMeans(object[,-column,drop=FALSE],na.rm=TRUE)
	y <- object[,column]-Ave
	x <- (object[,column]+Ave)/2

	plotWithHighlights(x,y,xlab=xlab,ylab=ylab,main=main,status=status, ...)
}


mdplot <- function(x,columns=c(1,2),xlab="Mean",ylab="Difference",main=NULL,...)
#	Mean-difference plot of two columns of a matrix
#	Gordon Smyth
#	16 March 2005. Last modified 13 April 2015.
{
	x <- as.matrix(x)
	columns <- as.integer(columns)[1:2]
	d <- x[,columns[2]]-x[,columns[1]]
	m <- (x[,columns[1]]+x[,columns[2]])/2
	if(is.null(main)) {
		cn <- colnames(x)[columns]
		if(is.null(cn)) cn <- paste("Column",columns)
		main <- paste(cn[2],"vs",cn[1])
	}
	plotWithHighlights(x=m,y=d,xlab=xlab,ylab=ylab,main=main,...)
}
