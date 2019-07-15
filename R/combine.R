# CBIND, RBIND AND MERGE

cbind.RGList <- function(..., deparse.level=1)
#  Combine RGList objects assuming same genelists
#  Gordon Smyth
#  27 June 2003. Last modified 6 Nov 2005.
{
	objects <- list(...)
	nobjects <- length(objects)
	out <- objects[[1]]
	other <- names(objects[[1]]$other)
	if(nobjects > 1)
	for (i in 2:nobjects) {
		out$R <- cbind(out$R,objects[[i]]$R)
		out$G <- cbind(out$G,objects[[i]]$G)
		out$Rb <- cbind(out$Rb,objects[[i]]$Rb)
		out$Gb <- cbind(out$Gb,objects[[i]]$Gb)
		out$weights <- cbind(out$weights,objects[[i]]$weights)
		out$targets <- rbind(out$targets,objects[[i]]$targets)
		for (a in other) out$other[[a]] <- cbind(out$other[[a]],objects[[i]]$other[[a]])
	}
	out
}

cbind.MAList <- function(..., deparse.level=1)
#  Combine MAList objects assuming same genelists
#  Gordon Smyth
#  27 June 2003. Last modified 6 Nov 2005.
{
	objects <- list(...)
	nobjects <- length(objects)
	out <- objects[[1]]
	other <- names(objects[[1]]$other)
	if(nobjects > 1)
	for (i in 2:nobjects) {
		out$M <- cbind(out$M,objects[[i]]$M)
		out$A <- cbind(out$A,objects[[i]]$A)
		out$weights <- cbind(out$weights,objects[[i]]$weights)
		out$targets <- rbind(out$targets,objects[[i]]$targets)
		for (a in other) out$other[[a]] <- cbind(out$other[[a]],objects[[i]]$other[[a]])
	}
	out
}

cbind.EListRaw <- cbind.EList <- function(..., deparse.level=1)
#  Combine EList objects assuming same genelists
#  Gordon Smyth
#  23 March 2009.  Last modified 5 June 2013.
{
	objects <- list(...)
	nobjects <- length(objects)
	out <- objects[[1]]
	other <- names(objects[[1]]$other)
	if(nobjects > 1)
	for (i in 2:nobjects) {
		out$E <- cbind(out$E,objects[[i]]$E)
		out$Eb <- cbind(out$Eb,objects[[i]]$Eb)
		out$weights <- cbind(out$weights,objects[[i]]$weights)
		out$targets <- rbind(out$targets,objects[[i]]$targets)
		out$design <- rbind(out$design,objects[[i]]$design)
		for (a in other) out$other[[a]] <- cbind(out$other[[a]],objects[[i]]$other[[a]])
	}
	out
}

rbind.RGList <- function(..., deparse.level=1)
#  Combine RGList objects assuming same array lists
#  Gordon Smyth
#  6 Dec 2003. Last modified 6 Nov 2005.
{
	objects <- list(...)
	nobjects <- length(objects)
	out <- objects[[1]]
	other <- names(objects[[1]]$other)
	if(nobjects > 1)
	for (i in 2:nobjects) {
		out$R <- rbind(out$R,objects[[i]]$R)
		out$G <- rbind(out$G,objects[[i]]$G)
		out$Rb <- rbind(out$Rb,objects[[i]]$Rb)
		out$Gb <- rbind(out$Gb,objects[[i]]$Gb)
		out$weights <- rbind(out$weights,objects[[i]]$weights)
		out$genes <- rbind(out$genes,objects[[i]]$genes)
		for (a in other) out$other[[a]] <- rbind(out$other[[a]],objects[[i]]$other[[a]])
	}
	out
}

rbind.MAList <- function(..., deparse.level=1)
#  Combine MAList objects assuming same array lists
#  Gordon Smyth
#  7 Dec 2003. Last modified 6 Nov 2005.
{
	objects <- list(...)
	nobjects <- length(objects)
	out <- objects[[1]]
	other <- names(objects[[1]]$other)
	if(nobjects > 1)
	for (i in 2:nobjects) {
		out$M <- rbind(out$M,objects[[i]]$M)
		out$A <- rbind(out$A,objects[[i]]$A)
		out$weights <- rbind(out$weights,objects[[i]]$weights)
		out$genes <- rbind(out$genes,objects[[i]]$genes)
		for (a in other) out$other[[a]] <- rbind(out$other[[a]],objects[[i]]$other[[a]])
	}
	out
}

rbind.EListRaw <- rbind.EList <- function(..., deparse.level=1)
#  Combine EList objects assuming same array lists
#  Gordon Smyth
#  23 March 2009.  Last modified 26 October 2010.
{
	objects <- list(...)
	nobjects <- length(objects)
	out <- objects[[1]]
	other <- names(objects[[1]]$other)
	am <- function(x) if(is.null(x)) NULL else as.matrix(x)
	if(nobjects > 1)
	for (i in 2:nobjects) {
		out$E <- rbind(am(out$E),am(objects[[i]]$E))
		out$Eb <- rbind(am(out$Eb),am(objects[[i]]$Eb))
		out$weights <- rbind(am(out$weights),am(objects[[i]]$weights))
		out$genes <- rbind(out$genes,objects[[i]]$genes)
		for (a in other) out$other[[a]] <- rbind(am(out$other[[a]]),am(objects[[i]]$other[[a]]))
	}
	out
}

makeUnique <- function(x)
#  Add characters to the elements of a character vector to make all values unique
#  Gordon Smyth
#  10 April 2003.  Last modified 28 June 2016.
{
	x <- as.character(x)
	tab <- table(x)
	tab <- tab[tab>1]
	lentab <- length(tab)
	if(lentab > 0) {
		u <- names(tab)
		for (i in 1:lentab) {
			n <- tab[i]
			x[x==u[i]] <- paste0(x[x==u[i]],formatC(1:n,width=1+floor(log10(n)),flag="0"))
		}
	}
	x
}

merge.RGList <- function(x,y,...)
#  Merge RGList y into x aligning by row names
#  Gordon Smyth
#  11 April 2003.  Last modified 28 Oct 2005.
{
	if(!is(y,"RGList")) stop("both x and y must be RGList objects")
	genes1 <- rownames(x$R)
	if(is.null(genes1)) genes1 <- rownames(x$G)
	if(is.null(genes1)) genes1 <- x$genes$ID
	genes2 <- rownames(y$R)
	if(is.null(genes2)) genes2 <- rownames(y$G)
	if(is.null(genes2)) genes2 <- y$genes$ID
	if(is.null(genes1) || is.null(genes2)) stop("Need row names to align on") 

	fields1 <- names(x)
	fields2 <- names(y)
	if(!identical(fields1,fields2)) stop("The two RGLists have different components")

	ord2 <- match(makeUnique(genes1), makeUnique(genes2))
	cbind(x,y[ord2,])
}

merge.MAList <- function(x,y,...)
#  Merge MAList y into x aligning by row names
#  Gordon Smyth
#  7 May 2004.  Last modified 28 Oct 2005.
{
	if(!is(y,"MAList")) stop("both x and y must be MAList objects")
	genes1 <- rownames(x$M)
	if(is.null(genes1)) genes1 <- rownames(x$A)
	if(is.null(genes1)) genes1 <- x$genes$ID
	genes2 <- rownames(y$M)
	if(is.null(genes2)) genes2 <- rownames(y$A)
	if(is.null(genes2)) genes2 <- y$genes$ID
	if(is.null(genes1) || is.null(genes2)) stop("Need row names to align on") 

	fields1 <- names(x)
	fields2 <- names(y)
	if(!identical(fields1,fields2)) stop("The two MALists have different components")

	ord2 <- match(makeUnique(genes1), makeUnique(genes2))
	cbind(x,y[ord2,])
}

merge.EListRaw <- function(x,y,...)
#  Merge EListRaw y into x aligning by row names
#  Gordon Smyth
#  9 May 2013.  Last modified 9 May 2013.
{
	if(!is(y,"EListRaw")) stop("both x and y must be EListRaw objects")
	genes1 <- rownames(x$E)
	if(is.null(genes1)) genes1 <- row.names(x$genes)
	if(is.null(genes1)) genes1 <- x$genes$ID
	genes2 <- rownames(y$E)
	if(is.null(genes2)) genes2 <- row.names(y$genes)
	if(is.null(genes2)) genes2 <- y$genes$ID
	if(is.null(genes1) || is.null(genes2)) stop("Need row names to align on") 

	fields1 <- names(x)
	fields2 <- names(y)
	if(!identical(fields1,fields2)) stop("The two MALists have different components")

	ord2 <- match(makeUnique(genes1), makeUnique(genes2))
	cbind(x,y[ord2,])
}

merge.EList <- function(x,y,...)
#  Merge EList y into x aligning by row names
#  Gordon Smyth
#  9 May 2013.  Last modified 9 May 2013.
{
	if(!is(y,"EList")) stop("both x and y must be EList objects")
	genes1 <- rownames(x$E)
	if(is.null(genes1)) genes1 <- row.names(x$genes)
	if(is.null(genes1)) genes1 <- x$genes$ID
	genes2 <- rownames(y$E)
	if(is.null(genes2)) genes2 <- row.names(y$genes)
	if(is.null(genes2)) genes2 <- y$genes$ID
	if(is.null(genes1) || is.null(genes2)) stop("Need row names to align on") 

	fields1 <- names(x)
	fields2 <- names(y)
	if(!identical(fields1,fields2)) stop("The two MALists have different components")

	ord2 <- match(makeUnique(genes1), makeUnique(genes2))
	cbind(x,y[ord2,])
}

