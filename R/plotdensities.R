#  PLOT DENSITIES

plotDensities <- function(object,...)
UseMethod("plotDensities")

plotDensities.RGList <- function(object,log=TRUE,group=NULL,col=NULL,main="RG Densities",bc.method="subtract",...)
#	Plot empirical single-channel densities
#	Original version by Natalie Thorne, 9 September 2003
#	Modified by Gordon Smyth.  Last modified 13 March 2015.
{
	object <- backgroundCorrect(object,method=bc.method)
	narray <- ncol(object)
	E <- cbind(object$R,object$G)

#	Add one to prevent taking logs of zero
#	However this means that plotDensities.RGList and plotDensities.MAList will give slightly
#	different plots on the same data even when the RG and MA objects represent the same data
	if(log) E <- log2(E+1)

	col2 <- col
	if(is.null(group)) {
		group2 <- factor(rep(1:2,c(narray,narray)),labels=c("R","G"))
		if(is.null(col2)) col2 <- c("red","green")
	} else {
		group <- rep(group,narray)
		group2 <- c(group,group)
	}

	plotDensities(object=E,group=group2,col=col2,main=main,...)
}

plotDensities.MAList <- function(object,log=TRUE,group=NULL,col=NULL,main="RG Densities",...)
#	Plot empirical single-channel densities
#	Original version by Natalie Thorne, 9 September 2003
#	Modified by Gordon Smyth.  Last modified 1 March 2015.
{
	narray <- ncol(object)
	E <- cbind(object$A+object$M/2, object$A-object$M/2)
	if(!log) E <- 2^E

	col2 <- col
	if(is.null(group)) {
		group2 <- factor(rep(1:2,c(narray,narray)),labels=c("R","G"))
		if(is.null(col2)) col2 <- c("red","green")
	} else {
		group <- rep(group,narray)
		group2 <- c(group,group)
	}
	
	plotDensities(object=E,group=group2,col=col2,main=main,...)
}

plotDensities.EListRaw <- function(object,log=TRUE,bc.method="subtract",...)
#	Gordon Smyth.
#	Created 23 March 2009.  Last modified 13 March 2015.
{
	object <- backgroundCorrect(object,method=bc.method)
	E <- object$E
	if(log) E <- log2(E+1)
	plotDensities(object=E,...)
}

plotDensities.EList <- function(object,log=TRUE,...)
#	Gordon Smyth.
#	Created 23 March 2009.  Last modified 13 March 2015.
{
	E <- object$E
	if(!log) E <- 2^E
	plotDensities(object=E,...)
}

plotDensities.default <- function(object,group=NULL,col=NULL,main=NULL,legend="topleft",...)
#	Plot empirical single-channel densities
#	Gordon Smyth
#	18 Nov 2013.  Last modified 13 March 2015.
{
#	Coerce object to matrix
	E <- as.matrix(object)
	narray <- ncol(E)

#	Check group
	if(is.null(group))  group <- colnames(E)
	if(is.null(group))  group <- 1:narray
	group <- as.factor(group)
	ngroup <- nlevels(group)

#	Check col
	if(is.null(col)) col <- 1:ngroup
	col <- rep(col,length=ngroup)

#	Check legend
	if(is.logical(legend)) {
		legend.position <- "topleft"
	} else {
		legend.position <- as.character(legend)
		legend <- TRUE
	}
	legend.position <- match.arg(legend.position,c("bottomright","bottom","bottomleft","left","topleft","top","topright","right","center"))

#	Expand cols to number of arrays
	arraycol <- group
	levels(arraycol) <- col
	arraycol <- as.vector(arraycol)

	npoint <- 512
	X <- Y <- matrix(0,npoint,narray)
	for (a in 1:ncol(E)) {
		d <- density(E[,a],n=npoint,na.rm=TRUE,...)
		X[,a] <- d$x
		Y[,a] <- d$y
	}
	matplot(X,Y,xlab="Intensity",ylab="Density",main=main,type="l",col=arraycol,lwd=2,lty=1)
	if(legend && ngroup>1) legend(legend.position,lwd=2,legend=levels(group),col=col)
	invisible(list(X=X,Y=Y))
}
