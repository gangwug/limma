#  SUBSET DATA SETS

subsetListOfArrays <- function(object,i,j,IJ,IX,I,JX)
#	Subsetting for list-like data objects
#	Gordon Smyth
#	11 Dec 2013
{
#	object,IJ,IX,I,JX are required arguments

#	Remove components that are absent or scalars
	len <- vapply(object,length,0)
	I <- intersect(I,names(len)[len>1L])

	if(missing(i)) {
		IX <- I <- character(0)
		if(missing(j)) IJ <- character(0)
	} else {
		if(is.character(i)) {
			i <- match(i, rownames(object))
			if(any(is.na(i))) stop("Subscript not found in rownames")
		}
	}
	if(missing(j)) {
		JX <- character(0)
	} else {
		if(is.character(j)) {
			j <- match(j, colnames(object))
			if(any(is.na(j))) stop("Subscript not found in colnames")
		}
	}
	
	for(a in IJ) object[[a]] <- object[[a]][i,j,drop=FALSE]
	for(a in IX) object[[a]] <- object[[a]][i, ,drop=FALSE]
	for(a in I ) object[[a]] <- object[[a]][i]
	for(a in JX) object[[a]] <- object[[a]][j, ,drop=FALSE]
	
	object
}

assign("[.RGList",
function(object, i, j)
#  Subsetting for RGList objects
#  Gordon Smyth
#  29 June 2003.  Last modified 20 July 2015.
{
	if(nargs() != 3) stop("Two subscripts required",call.=FALSE)

#	Recognized components
	IJ <- c("R","G","Rb","Gb","weights")
	IX <- c("genes")
	JX <- c("targets")
	I  <- character(0)
	object <- subsetListOfArrays(object,i,j,IJ=IJ,IX=IX,I=I,JX=JX)

	oc <- names(object$other)
	if(!missing(i) || !missing(j)) for(a in oc) object$other[[a]] <- object$other[[a]][i,j,drop=FALSE]
	
	object
})

assign("[.MAList",
function(object, i, j) 
#  Subsetting for MAList objects
#  Gordon Smyth
#  29 June 2003.  Last modified 20 July 2015.
{
	if(nargs() != 3) stop("Two subscripts required",call.=FALSE)

#	Recognized components
	IJ <- c("M","A","weights")
	IX <- c("genes")
	JX <- c("targets","design")
	I  <- character(0)
	object <- subsetListOfArrays(object,i,j,IJ=IJ,IX=IX,I=I,JX=JX)

	if(!missing(j) && !is.null(object$design) && !is.fullrank(object$design)) warning("subsetted design matrix is singular",call.=FALSE)

	oc <- names(object$other)
	if(!missing(i) || !missing(j)) for(a in oc) object$other[[a]] <- object$other[[a]][i,j,drop=FALSE]
	
	object
})

assign("[.EList",
function(object, i, j)
#  Subsetting for EList objects
#  Gordon Smyth
#  23 February 2009.  Last modified 20 July 2015.
{
	if(nargs() != 3) stop("Two subscripts required",call.=FALSE)

#	Recognized components
	IJ <- c("E","Eb","weights")
	IX <- c("genes")
	JX <- c("targets","design")
	I  <- character(0)
	object <- subsetListOfArrays(object,i,j,IJ=IJ,IX=IX,I=I,JX=JX)

	oc <- names(object$other)
	if(!missing(i) || !missing(j)) for(a in oc) object$other[[a]] <- object$other[[a]][i,j,drop=FALSE]
	
	object
})

assign("[.EListRaw", get("[.EList"))

assign("[.MArrayLM",
function(object, i, j)
#  Subsetting for MArrayLM objects
#  Gordon Smyth
#  26 April 2005. Last modified 20 July 2015.
{
	if(nargs() != 3) stop("Two subscripts required",call.=FALSE)

#	Recognized components
	IJ <- c("coefficients","stdev.unscaled","t","p.value","lods","weights")
	IX <- "genes"
	JX <- character(0)
	I  <- c("Amean","sigma","df.residual","df.prior","df.total","s2.post","F","F.p.value")
	JJ <- "cov.coefficients"
	XJ <- "contrasts"
	J  <- "var.prior"

#	After subsetting by columns, a contrast component should always be present
#	so that output is equivalent to that from contrasts.fit()
	if(!missing(j) && is.null(object$contrasts) && !is.null(object$coefficients)) {
		object$contrasts <- diag(ncol(object$coefficients))
		cn <- colnames(object$coefficients)
		dimnames(object$contrasts) <- list(Coefficient=cn,Contrast=cn)
	}

#	Ensure matrix or data.frame objects not dropped to vectors
	for (a in c(IJ,JJ,XJ)) if(!is.null(object[[a]])) object[[a]] <- as.matrix(object[[a]])
	for (a in c("targets","genes")) if(!is.null(object[[a]]) && is.null(dim(object[[a]]))) object[[a]] <- data.frame(object[[a]])

	object <- subsetListOfArrays(object,i,j,IJ=IJ,IX=IX,I=I,JX=JX)

#	Special treatment for JJ,XJ,J
	if(!missing(j)) {
		object$cov.coefficients <- object$cov.coefficients[j,j,drop=FALSE]
		object$contrasts <- object$contrasts[,j,drop=FALSE]
		object$var.prior <- object$var.prior[j]
	}

#	If columns have been subsetted, need to re-generate F
	if(!is.null(object[["F"]]) && !missing(j)) {
		F.stat <- classifyTestsF(object,fstat.only=TRUE)
		object$F <- as.vector(F.stat)
		df1 <- attr(F.stat,"df1")
		df2 <- attr(F.stat,"df2")
		if (df2[1] > 1e6) 
			object$F.p.value <- pchisq(df1*object$F,df1,lower.tail=FALSE)
		else
			object$F.p.value <- pf(object$F,df1,df2,lower.tail=FALSE)
	}

	object
})

