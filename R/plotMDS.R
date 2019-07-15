##  PLOTMDS.R

#	Class to hold multidimensional scaling output
setClass("MDS",representation("list"))

setMethod("show","MDS",function(object) {
	cat("An object of class MDS\n")
	print(unclass(object))
})

plotMDS <- function(x,...) UseMethod("plotMDS")

plotMDS.MDS <- function(x,labels=NULL,pch=NULL,cex=1,dim.plot=NULL,xlab=NULL,ylab=NULL,...)
#	Method for MDS objects
#	Create a new plot using MDS coordinates or distances previously created
#	Gordon Smyth and Yifang Hu
#	21 May 2011.  Last modified 29 December 2014
{
#	Check labels
	if(is.null(labels) & is.null(pch)) {
		labels <- colnames(x$distance.matrix)
		if(is.null(labels)) labels <- 1:length(x$x)
	}

#	Are new dimensions requested?
	if(is.null(dim.plot)) {
		dim.plot <- x$dim.plot
	} else {
		if(any(dim.plot != x$dim.plot)) {
			ndim <- max(dim.plot)
			if(ndim > ncol(x$cmdscale.out)) x$cmdscale.out <- cmdscale(as.dist(x$distance.matrix),k=ndim)
			x$x <- x$cmdscale.out[,dim.plot[1]]
			x$y <- x$cmdscale.out[,dim.plot[2]]
		}
	}

#	Axis labels
	if(is.null(x$axislabel)) x$axislabel <- "Principal Coordinate"
	if(is.null(xlab)) xlab <- paste(x$axislabel,dim.plot[1])
	if(is.null(ylab)) ylab <- paste(x$axislabel,dim.plot[2])

#	Make the plot
	if(is.null(labels)){
#		Plot symbols instead of text
		plot(x$x, x$y, pch = pch, xlab = xlab, ylab = ylab, cex = cex, ...)
	} else {
#		Plot text.
		labels <- as.character(labels)
#		Need to estimate width of labels in plot coordinates.
#		Estimate will be ok for default plot width, but maybe too small for smaller plots.
		StringRadius <- 0.01*cex*nchar(labels)
		left.x <- x$x-StringRadius
		right.x <- x$x+StringRadius
		plot(c(left.x, right.x), c(x$y, x$y), type = "n", xlab = xlab, ylab = ylab, ...)
		text(x$x, x$y, labels = labels, cex = cex, ...)
	}

	invisible(x)
}

plotMDS.default <- function(x,top=500,labels=NULL,pch=NULL,cex=1,dim.plot=c(1,2),ndim=max(dim.plot),gene.selection="pairwise",xlab=NULL,ylab=NULL,plot=TRUE,...)
#	Multi-dimensional scaling with top-distance
#	Di Wu and Gordon Smyth
#	19 March 2009.  Last modified 6 Oct 2016
{
#	Check x
	x <- as.matrix(x)
	nsamples <- ncol(x)
	if(nsamples < 3) stop(paste("Only",nsamples,"columns of data: need at least 3"))
	cn <- colnames(x)
#	Remove rows with missing or Inf values
	bad <- rowSums(is.finite(x)) < nsamples
	if(any(bad)) x <- x[!bad,,drop=FALSE]
	nprobes <- nrow(x)

#	Check top
	top <- min(top,nprobes)

#	Check labels and pch
	if(is.null(pch) & is.null(labels)) {
		labels <- colnames(x)
		if(is.null(labels)) labels <- 1:nsamples
	}
	if(!is.null(labels)) labels <- as.character(labels)

#	Check dim.plot
	dim.plot <- unique(as.integer(dim.plot))
	if(length(dim.plot) != 2L) stop("dim.plot must specify two dimensions to plot")

#	Check dim
	if(ndim < 2L) stop("Need at least two dim.plot")
	if(nsamples < ndim) stop("ndim is greater than number of samples")
	if(nprobes < ndim) stop("ndim is greater than number of rows of data")

#	Check gene.selection
	gene.selection <- match.arg(gene.selection,c("pairwise","common"))

#	Distance matrix from pairwise leading fold changes
	dd <- matrix(0,nrow=nsamples,ncol=nsamples,dimnames=list(cn,cn))
	if(gene.selection=="pairwise") {
#		Distance measure is mean of top squared deviations for each pair of arrays
		topindex <- nprobes-top+1L
		for (i in 2L:(nsamples))
		for (j in 1L:(i-1L))
			dd[i,j]=sqrt(mean(sort.int((x[,i]-x[,j])^2,partial=topindex)[topindex:nprobes]))
		axislabel <- "Leading logFC dim"
	} else {
#		Same genes used for all comparisons
		if(nprobes > top) {
			s <- rowMeans((x-rowMeans(x))^2)
			o <- order(s,decreasing=TRUE)
			x <- x[o[1L:top],,drop=FALSE]
		}
		for (i in 2L:(nsamples))
			dd[i,1L:(i-1L)]=sqrt(colMeans((x[,i]-x[,1:(i-1),drop=FALSE])^2))
		axislabel <- "Principal Component"
	}

#	Multi-dimensional scaling
	a1 <- suppressWarnings(cmdscale(as.dist(dd),k=ndim))

#	Make MDS object and call plotMDS method
	mds <- new("MDS",list(dim.plot=dim.plot,distance.matrix=dd,cmdscale.out=a1,top=top,gene.selection=gene.selection))
	if(dim.plot[1] > ncol(a1)) {
		mds$x <- rep.int(0,nsamples)
		warning(paste("dimension",dim.plot[1],"is degenerate or all zero"))
	} else
		mds$x <- a1[,dim.plot[1]]
	if(dim.plot[2] > ncol(a1)) {
		mds$y <- rep.int(0,nsamples)
		warning(paste("dimension",dim.plot[2],"is degenerate or all zero"))
	} else
		mds$y <- a1[,dim.plot[2]]
	mds$top <- top
	mds$axislabel <- axislabel
	if(plot)
		plotMDS(mds,labels=labels,pch=pch,cex=cex,xlab=xlab,ylab=ylab,...)
	else
		mds
}
