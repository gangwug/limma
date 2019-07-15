#  avearrays.R

avearrays <- function(x,ID=NULL,weights=NULL) UseMethod("avearrays")
#  24 Sept 2010

avearrays.default <- function(x,ID=colnames(x),weights=NULL)
#	Average over technical replicate columns, for matrices
#	Gordon Smyth
#	Created 24 Sept 2010. Last modified 13 Aug 2015.
{
	if(is.null(x)) return(NULL)
	if(is.vector(x)) stop("x must be a matrix")
	x <- as.matrix(x)
	if(is.null(ID)) stop("No sample IDs")
	ID <- as.character(ID)
	if(mode(x)=="character") {
		d <- duplicated(ID)
		if(!any(d)) return(x)
		y <- x[,!d,drop=FALSE]
		return(y)
	}
	ID <- factor(ID,levels=unique(ID))
	if(is.null(weights)) {
		y <- t(rowsum(t(x),ID,reorder=FALSE,na.rm=TRUE))
		n <- t(rowsum(t(1L-is.na(x)),ID,reorder=FALSE,na.rm=TRUE))
		y <- y/n
	} else {
		design <- model.matrix(~0+ID)
		y <- lmFit(x,design,weights=weights)$coefficients
	}
	y
}

avearrays.MAList <- function(x,ID=colnames(x),weights=x$weights)
#	Average over technical replicate columns for MAList objects
#	Gordon Smyth
#	24 Sept 2010.  Last modified 24 Sep 2010.
{
	ID <- as.character(ID)
	d <- duplicated(ID)
	if(!any(d)) return(x)
	y <- x[,!d]
	y$M <- avearrays(x$M,ID,weights=weights)
	y$A <- avearrays(x$A,ID,weights=weights)
	y$weights <- avearrays(x$weights,ID,weights=weights)
	other <- names(x$other)
	for (a in other) y$other[[a]] <- avearrays(x$other[[a]],ID=ID,weights=weights)
	y
}

avearrays.EList <- function(x,ID=colnames(x),weights=x$weights)
#	Average over irregular replicate columns for EList objects
#	Gordon Smyth
#	24 Sept 2010.  Last modified 27 Oct 2010.
{
	d <- duplicated(ID)
	if(!any(d)) return(x)
	y <- x[,!d]
	y$E <- avearrays(x$E,ID,weights=weights)
	y$weights <- avearrays(x$weights,ID,weights=weights)
	other <- names(x$other)
	for (a in other) y$other[[a]] <- avearrays(x$other[[a]],ID=ID,weights=weights)
	y
}
