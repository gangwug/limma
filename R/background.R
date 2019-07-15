#  BACKGROUND.R

#  BACKGROUND CORRECTION

backgroundCorrect <- function(RG, method="auto", offset=0, printer=RG$printer, normexp.method="saddle", verbose=TRUE)
{
#	Apply background correction to microarray data
#	Gordon Smyth
#	12 April 2003.  Last modified 30 August 2010.

	if(is.null(method)) method <- "auto"

	if(is(RG,"EListRaw")) {
		RG$E <- backgroundCorrect.matrix(RG$E,RG$Eb,method=method,offset=offset,normexp.method=normexp.method,verbose=verbose)
		RG$Eb <- NULL
		return(RG)
	}

	if(class(RG)=="list" && !is.null(RG$R)) RG <- new("RGList",RG)

	if(is(RG,"RGList")) {
		RG$R <- backgroundCorrect.matrix(RG$R,RG$Rb,method=method,offset=offset,normexp.method=normexp.method,verbose=verbose)
		RG$Rb <- NULL
		RG$G <- backgroundCorrect.matrix(RG$G,RG$Gb,method=method,offset=offset,normexp.method=normexp.method,verbose=verbose)
		RG$Gb <- NULL
		return(RG)
	}

	RG <- as.matrix(RG)
	if(!(method %in% c("none","normexp"))) stop(method,"correction requires background intensities")
	RG <- backgroundCorrect.matrix(RG,method=method,offset=offset,normexp.method=normexp.method,verbose=verbose)
	return(RG)
}


backgroundCorrect.matrix <- function(E, Eb=NULL, method="auto", offset=0, printer=NULL, normexp.method="saddle", verbose=TRUE)
{
#	Apply background correction to microarray data
#	in the form of foreground and background matrices
#	Gordon Smyth
#	30 August 2010.  Last modified 30 August 2010.

	E <- as.matrix(E)
	if(method=="auto") method <- ifelse(is.null(Eb),"normexp","subtract")
	method <- match.arg(method, c("none","subtract","half","minimum","movingmin","edwards","normexp"))

	if(is.null(Eb)) {
		if(method %in% c("subtract","half","minimum","movingmin","edwards")) method <- "none"
	} else {
		Eb <- as.matrix(Eb)
	}

	switch(method,
	subtract={
		E <- E-Eb
	},
	half={
		E <- pmax(E-Eb, 0.5)
	},
	minimum={
		E <- E-Eb
		for (slide in 1:ncol(E)) {
			i <- E[,slide] < 1e-18
			if(any(i,na.rm=TRUE)) {
				m <- min(E[!i,slide],na.rm=TRUE)
				E[i,slide] <- m/2
			}
		}
	},
	movingmin={
		if(is.null(printer)) {
			E <- E-ma3x3.matrix(Eb,FUN=min,na.rm=TRUE)
		} else {
			E <- E-ma3x3.spottedarray(Eb,printer=printer,FUN=min,na.rm=TRUE)
		}
	},
	edwards={
#		Log-linear interpolation for dull spots as in Edwards (2003).
#		The threshold values (delta) are chosen such that the number of
#		spots with (0 < R-Rb < delta) is f=10% of the number spots
#		with (R-Rb <= 0) for each channel and array.
#		Note slight change from Edwards (2003).
		one <- matrix(1,nrow(E),1)
		delta.vec <- function(d, f=0.1) {
			quantile(d, mean(d<1e-16,na.rm=TRUE)*(1+f), na.rm=TRUE)
		}
		sub <- E-Eb
		delta <- one %*% apply(sub, 2, delta.vec)
		E <- ifelse(sub < delta, delta*exp(1-(Eb+delta)/E), sub)
	},
	normexp={
		if(!is.null(Eb)) E <- E-Eb
		for (j in 1:ncol(E)) {
			if(verbose) cat("Array",j)
			x <- E[,j]
			out <- normexp.fit(x,method=normexp.method)
			E[,j] <- normexp.signal(out$par,x)
			if(verbose) cat(" corrected\n")
		}
	}
	)
	if(offset) E <- E+offset
	E
}

ma3x3.matrix <- function(x,FUN=mean,na.rm=TRUE,...)
#	2-dimensional moving average for 3x3 blocks
#	Gordon Smyth
#	11 April 2004
{
#	Pad out x with NA so that original values have 8 neighbors
	d1 <- nrow(x)
	d2 <- ncol(x)
	y <- matrix(NA,d1+2,d2+2)
	y[1+(1:d1),1+(1:d2)] <- x

#	Index vector for original values
	i <- 1:length(y)
	dim(i) <- dim(y)
	i <- i[1+(1:d1),1+(1:d2)]
	dim(i) <- NULL

#	Rows are original obs, columns are neighbors
	x <- matrix(x,d1*d2,9)
	ry <- nrow(y)
	x[,1] <- y[i-ry-1]
	x[,2] <- y[i-ry]
	x[,3] <- y[i-ry+1]
	x[,4] <- y[i-1]
	x[,6] <- y[i+1]
	x[,7] <- y[i+ry-1]
	x[,8] <- y[i+ry]
	x[,9] <- y[i+ry+1]

	y <- apply(x,MARGIN=1,FUN=FUN,na.rm=na.rm,...)
	dim(y) <- c(d1,d2)
	y
}

ma3x3.spottedarray <- function(x,printer,FUN=mean,na.rm=TRUE,...)
#	Gordon Smyth
#	11 April 2004
{
	x <- as.matrix(x)
	narrays <- ncol(x)
	gr <- printer$ngrid.r
	gc <- printer$ngrid.c
	sr <- printer$nspot.r
	sc <- printer$nspot.c
	dim(x) <- c(sc, sr, gc, gr, narrays)
	x <- aperm(x, perm = c(2, 4, 1, 3, 5))
	dim(x) <- c(gr * sr, gc * sc, narrays)
	for (j in 1:narrays) x[,,j] <- ma3x3.matrix(x[,,j],FUN=FUN,na.rm=TRUE,...)
	dim(x) <- c(sr, gr, sc, gc, narrays)
	x <- aperm(x, perm = c(3, 1, 4, 2, 5))
	dim(x) <- c(sc*sr*gc*gr, narrays)
	x
}

