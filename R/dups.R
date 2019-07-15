#  DUPS.R
#  Functions to handle duplicate spots or blocking

unwrapdups <- function(M,ndups=2,spacing=1) {
#	Unwrap M matrix for a series of experiments so that all spots for a given gene are in one row
#	Gordon Smyth
#	18 Jan 2002. Last revised 2 Nov 2002.

	if(ndups==1) return(M)
	M <- as.matrix(M)
	nspots <- dim(M)[1]
	nslides <- dim(M)[2]
	ngroups <- nspots / ndups / spacing
	dim(M) <- c(spacing,ndups,ngroups,nslides)
	M <- aperm(M,perm=c(1,3,2,4))
	dim(M) <- c(spacing*ngroups,ndups*nslides)
	M
}

uniquegenelist <- function(genelist,ndups=2,spacing=1) {
#	Eliminate entries in genelist for duplicate spots
#	Gordon Smyth
#	2 Nov 2002.  Last revised 10 Jan 2005

	if(ndups <= 1) return(genelist)
	i <- drop(unwrapdups(1:NROW(genelist),ndups=ndups,spacing=spacing)[,1])
	if(is.null(dim(genelist)))
		return(genelist[i])
	else
		return(genelist[i,,drop=FALSE])
}

duplicateCorrelation <- function(object,design=NULL,ndups=2,spacing=1,block=NULL,trim=0.15,weights=NULL)
#	Estimate the correlation between duplicates given a series of arrays
#	Gordon Smyth
#	25 Apr 2002. Last revised 17 Aug 2016.
{
#	Extract components from y
	y <- getEAWP(object)
	M <- y$exprs
	ngenes <- nrow(M)
	narrays <- ncol(M)

#	Check design matrix
	if(is.null(design)) design <- y$design
	if(is.null(design))
		design <- matrix(1,ncol(y$exprs),1)
	else {
		design <- as.matrix(design)
		if(mode(design) != "numeric") stop("design must be a numeric matrix")
	}
	if(nrow(design) != narrays) stop("Number of rows of design matrix does not match number of arrays")
	ne <- nonEstimable(design)
	if(!is.null(ne)) cat("Coefficients not estimable:",paste(ne,collapse=" "),"\n")
	nbeta <- ncol(design)

#	Weights and spacing arguments can be specified in call or stored in y
#	Precedence for these arguments is
#	1. Specified in function call
#	2. Stored in object
#	3. Default values
	if(missing(ndups) && !is.null(y$printer$ndups)) ndups <- y$printer$ndups
	if(missing(spacing) && !is.null(y$printer$spacing)) spacing <- y$printer$spacing
	if(missing(weights) && !is.null(y$weights)) weights <- y$weights

#	Check weights
	if(!is.null(weights)) {
		weights <- asMatrixWeights(weights,dim(M))
		weights[weights <= 0] <- NA
		M[!is.finite(weights)] <- NA
	}

#	Setup spacing or blocking arguments
	if(is.null(block)) {
		if(ndups<2) {
			warning("No duplicates: correlation between duplicates not estimable")
			return( list(cor=NA,cor.genes=rep(NA,nrow(M))) )
		}
		if(is.character(spacing)) {
			if(spacing=="columns") spacing <- 1
			if(spacing=="rows") spacing <- object$printer$nspot.c
			if(spacing=="topbottom") spacing <- nrow(M)/2
		}
		Array <- rep(1:narrays,rep(ndups,narrays))
	} else {
		ndups <- 1
		nspacing <- 1
		Array <- block
	}

#	Unwrap data to get all data for a gene in one row
	if(is.null(block)) {
		M <- unwrapdups(M,ndups=ndups,spacing=spacing)
		ngenes <- nrow(M)
		if(!is.null(weights)) weights <- unwrapdups(weights,ndups=ndups,spacing=spacing)
		design <- design %x% rep(1,ndups)
	}

	if(!requireNamespace("statmod",quietly=TRUE)) stop("statmod package required but is not installed")
	rho <- rep(NA,ngenes)
	nafun <- function(e) NA
	for (i in 1:ngenes) {
		y <- drop(M[i,])
		o <- is.finite(y)
		A <- factor(Array[o])
		nobs <- sum(o)
		nblocks <- length(levels(A))
		if(nobs>(nbeta+2) && nblocks>1 && nblocks<nobs-1) {
			y <- y[o]
			X <- design[o,,drop=FALSE]
			Z <- model.matrix(~0+A)
			if(!is.null(weights)) {
				w <- drop(weights[i,])[o]
				s <- tryCatch(statmod::mixedModel2Fit(y,X,Z,w,only.varcomp=TRUE,maxit=20)$varcomp,error=nafun)
			} else
				s <- tryCatch(statmod::mixedModel2Fit(y,X,Z,only.varcomp=TRUE,maxit=20)$varcomp,error=nafun)
			if(!is.na(s[1])) rho[i] <- s[2]/sum(s)
		}
	}
	arho <- atanh(pmax(-1,rho))
	mrho <- tanh(mean(arho,trim=trim,na.rm=TRUE))
	list(consensus.correlation=mrho,cor=mrho,atanh.correlations=arho)
}

avedups <- function(x,ndups,spacing,weights) UseMethod("avedups")

avedups.default <- function(x,ndups=2,spacing=1,weights=NULL)
#	Average over duplicate spots, for matrices or vectors
#	Gordon Smyth
#	6 Apr 2006.
{
	if(ndups==1) return(x)
	if(is.null(x)) return(NULL)
	x <- as.matrix(x)
	nspots <- dim(x)[1]
	nslides <- dim(x)[2]
	rn <- rownames(x)
	cn <- colnames(x)
	ngroups <- nspots / ndups / spacing
	dim(x) <- c(spacing,ndups,ngroups*nslides)
	x <- aperm(x,perm=c(2,1,3))
	if(mode(x)=="character")
		x <- x[1,,]
	else {
		if(is.null(weights))
			x <- colMeans(x,na.rm=TRUE)
		else {
			weights <- as.matrix(weights)
			dim(weights) <- c(spacing,ndups,ngroups*nslides)
			weights <- aperm(weights,perm=c(2,1,3))
			weights[is.na(weights) | is.na(x)] <- 0
			weights[weights<0] <- 0
			x <- colSums(weights*x,na.rm=TRUE)/colSums(weights)
		}
	}
	dim(x) <- c(spacing*ngroups,nslides)
	colnames(x) <- cn
	rownames(x) <- avedups(rn,ndups=ndups,spacing=spacing)
	x
}

avedups.MAList <- function(x,ndups=x$printer$ndups,spacing=x$printer$spacing,weights=x$weights)
#	Average over duplicate spots for MAList objects
#	Gordon Smyth
#	6 Apr 2006.
{
	if(is.null(ndups) || is.null(spacing)) stop("Must specify ndups and spacing")
	y <- x
	y$M <- avedups(x$M,ndups=ndups,spacing=spacing,weights=weights)
	y$A <- avedups(x$A,ndups=ndups,spacing=spacing,weights=weights)
	other <- names(x$other)
	for (a in other) object$other[[a]] <- avedups(object$other[[a]],ndups=ndups,spacing=spacing,weights=weights)
	y$weights <- avedups(x$weights,ndups=ndups,spacing=spacing)
	y$genes <- uniquegenelist(x$genes,ndups=ndups,spacing=spacing)
	y$printer <- NULL
	y
}

avedups.EList <- function(x,ndups=x$printer$ndups,spacing=x$printer$spacing,weights=x$weights)
#	Average over duplicate spots for EList objects
#	Gordon Smyth
#	2 Apr 2010.
{
	if(is.null(ndups) || is.null(spacing)) stop("Must specify ndups and spacing")
	y <- x
	y$E <- avedups(x$E,ndups=ndups,spacing=spacing,weights=weights)
	other <- names(x$other)
	for (a in other) object$other[[a]] <- avedups(object$other[[a]],ndups=ndups,spacing=spacing,weights=weights)
	y$weights <- avedups(x$weights,ndups=ndups,spacing=spacing)
	y$genes <- uniquegenelist(x$genes,ndups=ndups,spacing=spacing)
	y$printer <- NULL
	y
}

avereps <- function(x,...)
#	4 June 2008
	UseMethod("avereps")

avereps.default <- function(x,ID=rownames(x),...)
#	Average over irregular replicate spots, for matrices or vectors
#	Gordon Smyth
#	Created 3 June 2008.  Last modified 1 Dec 2010.
#	Revised 19 Aug 2009 following suggestions from Axel Klenk.
#	Revised 28 March 2010 following suggestion from Michael Lawrence.
{
	if(is.null(x)) return(NULL)
	x <- as.matrix(x)
	if(is.null(ID)) stop("No probe IDs")
	ID <- as.character(ID)
	if(mode(x)=="character") {
		d <- duplicated(ID)
		if(!any(d)) return(x)
		y <- x[!d,,drop=FALSE]
		return(y)
	}
	ID <- factor(ID,levels=unique(ID))
#	rowsum(x,ID,reorder=FALSE,na.rm=TRUE)/as.vector(table(ID))
	y <- rowsum(x,ID,reorder=FALSE,na.rm=TRUE)
	n <- rowsum(1L-is.na(x),ID,reorder=FALSE)
	y/n
}

avereps.MAList <- function(x,ID=NULL,...)
#	Average over irregular replicate spots for MAList objects
#	Gordon Smyth
#	3 June 2008.  Last modified 8 Sep 2010.
{
	if(is.null(ID)) {
		ID <- x$genes$ID
		if(is.null(ID)) ID <- rownames(x)
		if(is.null(ID)) stop("Cannot find probe IDs")
	}
	y <- x
	y$M <- avereps(x$M,ID=ID)
	y$A <- avereps(x$A,ID=ID)
	other <- names(x$other)
	for (a in other) y$other[[a]] <- avereps(x$other[[a]],ID=ID)
	y$weights <- avereps(x$weights,ID=ID)
	y$genes <- x$genes[!duplicated(ID),]
	y$printer <- NULL
	y
}

avereps.EList <- function(x,ID=NULL,...)
#	Average over irregular replicate probes for EList objects
#	Gordon Smyth
#	2 April 2010.  Last modified 20 May 2011.
{
	if(is.null(ID)) {
		ID <- x$genes$ID
		if(is.null(ID)) ID <- rownames(x)
		if(is.null(ID)) stop("Cannot find probe IDs")
	}
	y <- x
	y$E <- avereps(x$E,ID=ID)
	other <- names(x$other)
	for (a in other) y$other[[a]] <- avereps(x$other[[a]],ID=ID)
	y$weights <- avereps(x$weights,ID=ID)
	y$genes <- x$genes[!duplicated(ID),]
	y$printer <- NULL
	y
}

avereps.RGList <- function(x,ID=NULL,...)
#	Warn users that averaging should not be applied prior to normalization
#	Gordon Smyth
#	2 December 2013.
{
	stop("avereps should not be applied to an RGList object")
	invisible()
}


avereps.EListRaw <- function(x,ID=NULL,...)
#	Warn users that averaging should not be applied prior to normalization
#	Gordon Smyth
#	2 December 2013.
{
	stop("avereps should not be applied to an EListRaw object")
	invisible()
}
