normexp.fit.detection.p <- function(x,detection.p="Detection")
#  Estimate normexp parameters using negative control probes which are derived from probes' detection p values
#  Wei Shi and Gordon Smyth
#  Created 27 October 2010.  Modified 27 May 2013.
{
	if(is(x,"EListRaw")){
		if(is.character(detection.p)){
			other.colnames <- tolower(names(x[["other"]]))
			detection.index <- which(other.colnames %in% tolower(detection.p))
			if(length(detection.index)!=1)
				stop("Detection p values not found in the data.")
			detection.p <- x[["other"]][[detection.index]]
		} else {
			detection.p <- as.matrix(detection.p)
		}
		x <- as.matrix(x$E)
	} else {
		if(is.character(detection.p))
			stop("detection.p must be a numeric matrix (unless x is an EListRaw)")
		x <- as.matrix(x)
		detection.p <- as.matrix(detection.p)
	}
	
	if(!all(dim(x) == dim(detection.p)))
		stop("The supplied detection p value data do not have the same dimension as that of the intensity data.")
	
	narrays <- ncol(x)

#	Check whether pvalues are actually 1-pvalues
	y <- x[,1]
	p <- detection.p[,1]	
	if(p[which.max(y)] < p[which.min(y)])
		detection.p <- 1-detection.p

	mu <- sigma <- rep(NA, narrays)
	for(i in 1:narrays){
		y <- x[,i]
		p <- detection.p[,i]
		o <- order(p,y)
		y <- y[o]
		p <- p[o]
		j <- which(!duplicated(p))[-1]
		ync <- (y[j]+y[j-1])/2
		d <- p[j]-p[j-1]
#		if(any(d<0)) stop("detection p-values are not monotonic in the expression values for array",i)
		freq <- d/min(d)
		n <- sum(freq)
		mu[i] <- weighted.mean(ync,freq)
		v <- (ync-mu[i])^2
		sigma[i] <- sqrt(weighted.mean(v,freq)*n/(n-1))
	}
	alpha <- pmax(colMeans(x,na.rm=TRUE)-mu,10)
	cbind(mu=mu,logsigma=log(sigma),logalpha=log(alpha))
}

normexp.fit.control <- function(x,status=NULL,negctrl="negative",regular="regular",robust=FALSE)
#  Estimate normexp parameters using negative control probes
#  Wei Shi and Gordon Smyth
#  Created 17 April 2009. Last modified 14 January 2015.
{
	if(is(x, "EListRaw")) {
		if(is.null(status)) status <- x$genes$Status
		x <- x$E
	}
	x <- as.matrix(x)
	if(is.null(status)) stop("Probe status not found")

	xr <- x[tolower(status)==tolower(regular),,drop=FALSE]
	if(nrow(xr)==0) stop("No regular probes found")
	xn <- x[tolower(status)==tolower(negctrl),,drop=FALSE]
	if(nrow(xn)<2) stop("Fewer than two negative control probes found")

	if(robust) {
		if(!requireNamespace("MASS",quietly=TRUE)) stop("MASS package required but is not available")
		narrays <- ncol(xn)
		m <- s <- rep(0,narrays)
	#	Robustness is judged on the log-scale, assumed normal
		for (j in 1:ncol(xn)) {
			h <- MASS::huber(log(xn[,j]))
			m[j] <- h$mu
			s[j] <- h$s
		}
	#	Means and standard deviation are converted back to log-normal
		mu <- exp(m+s^2/2)
		omega <- exp(s^2)
		sigma <- sqrt(omega*(omega-1))*exp(m)
	} else {
		mu <- colMeans(xn,na.rm=TRUE)
		sigma <- sqrt(rowSums((t(xn)-mu)^2,na.rm=TRUE)/(nrow(xn)-1))
	}
	alpha <- pmax(colMeans(xr,na.rm=TRUE)-mu,10)
	cbind(mu=mu,logsigma=log(sigma),logalpha=log(alpha))
}

nec <- function(x,status=NULL,negctrl="negative",regular="regular",offset=16,robust=FALSE,detection.p="Detection")
#	Normexp background correction aided by negative controls.
#	Wei Shi and Gordon Smyth
#	Created 27 September 2010. Last modified 15 Nov 2015.
{
	if(is(x, "EListRaw")) {
		if(!is.null(x$Eb)) {
			x$E <- x$E-x$Eb
			x$Eb <- NULL
		}
		if(is.null(status)) status <- x$genes$Status
		if(any(tolower(status) %in% tolower(negctrl))) {
			normexp.par <- normexp.fit.control(x,status,negctrl,regular,robust)
		} else {
			normexp.par <- normexp.fit.detection.p(x,detection.p)
			message("Note: inferring mean and variance of negative control probe intensities from the detection p-values.")
		}
		for(i in 1:ncol(x)) x$E[,i] <- normexp.signal(normexp.par[i,], x$E[,i])
		x$E <- x$E + offset
	} else {
		x <- as.matrix(x)
		if(any(tolower(status) %in% tolower(negctrl))){
			normexp.par <- normexp.fit.control(x,status,negctrl,regular,robust)
		} else {
			normexp.par <- normexp.fit.detection.p(x,detection.p)
		}
		for(i in 1:ncol(x)) x[,i] <- normexp.signal(normexp.par[i,], x[,i])
		x <- x + offset
	}
	x
}

neqc <- function(x, status=NULL, negctrl="negative", regular="regular", offset=16, robust=FALSE, detection.p="Detection", ...)
#	Normexp background correction and quantile normalization using control probes
#	Wei Shi and Gordon Smyth
#	Created 17 April 2009. Last modified 27 October 2010.
{
	x.bg <- nec(x,status,negctrl,regular,offset,robust,detection.p)
	if(is(x.bg, "EListRaw")) {
		y <- normalizeBetweenArrays(x.bg, method="quantile", ...)
		if(is.null(status))
			status <- y$genes$Status
		if(!is.null(status)){
			y <- y[tolower(status) == tolower(regular), ]
			y$genes$Status <- NULL
		}
	} else {
		x.bg <- as.matrix(x.bg)
		y <- log2(normalizeBetweenArrays(x.bg, method="quantile", ...))
		if(!is.null(status))
			y <- y[tolower(status) == tolower(regular), ]
	}
	y
}


