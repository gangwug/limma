#  WITHIN ARRAY NORMALIZATION

MA.RG <- function(object, bc.method="subtract", offset=0) {
#	Convert RGList to MAList
#	Gordon Smyth
#	2 March 2003.  Last revised 9 Dec 2004.

	if(is.null(object$R) || is.null(object$G)) stop("Object doesn't contain R and G components")
	object <- backgroundCorrect(object, method=bc.method, offset=offset)
	R <- object$R
	G <- object$G

#	Log
	R[R <= 0] <- NA
	G[G <= 0] <- NA
	R <- log(R,2)
	G <- log(G,2)
	
#	Minus and Add
	object$R <- object$G <- object$Rb <- object$Gb <- object$other <- NULL
	object$M <- as.matrix(R-G)
	object$A <- as.matrix((R+G)/2)
	new("MAList",unclass(object))
}

RG.MA <- function(object) {
#	Convert MAList to RGList
#	Gordon Smyth
#	5 September 2003.  Last modified 9 Dec 2004.

	object$R <- 2^( object$A+object$M/2 )
	object$G <- 2^( object$A-object$M/2 )
	object$M <- NULL
	object$A <- NULL
	new("RGList",unclass(object))
}

normalizeWithinArrays <- function(object,layout=object$printer,method="printtiploess",weights=object$weights,span=0.3,iterations=4,controlspots=NULL,df=5,robust="M",bc.method="subtract",offset=0)
#	Within array normalization
#	Gordon Smyth
#	2 March 2003.  Last revised 16 March 2013.
{
#	Check input arguments
#	and get non-intensity dependent methods out of the way
	if(!is(object,"MAList")) object <- MA.RG(object,bc.method=bc.method,offset=offset)
	choices <- c("none","median","loess","printtiploess","composite","control","robustspline")
	method <- match.arg(method,choices)
	if(method=="none") return(object)
	if(is.vector(object$M)) object$M <- as.matrix(object$M)
	nprobes <- nrow(object$M)
	narrays <- ncol(object$M)
	if(!is.null(weights)) weights <- asMatrixWeights(weights,dim=c(nprobes,narrays))
	if(method=="median") {
		if(is.null(weights))
			for (j in 1:narrays) object$M[,j] <- object$M[,j] - median(object$M[,j],na.rm=TRUE)
		else
			for (j in 1:narrays) object$M[,j] <- object$M[,j] - weighted.median(object$M[,j],weights[,j],na.rm=TRUE)
		return(object)
	}

#	All remaining methods use regression of M-values on A-values
	if(is.vector(object$A)) object$A <- as.matrix(object$A)
	if(nrow(object$A) != nprobes) stop("row dimension of A doesn't match M")
	if(ncol(object$A) != narrays) stop("col dimension of A doesn't match M")
	switch(method,
		loess = {
			for (j in 1:narrays) {
				y <- object$M[,j]
				x <- object$A[,j]
				w <- weights[,j]
				object$M[,j] <- loessFit(y,x,w,span=span,iterations=iterations)$residuals
			}
		},
		printtiploess = {
			if(is.null(layout)) stop("Layout argument not specified")
			ngr <- layout$ngrid.r
			ngc <- layout$ngrid.c
			nspots <- layout$nspot.r * layout$nspot.c
			nprobes2 <- ngr*ngc*nspots
			if(nprobes2 != nprobes) stop("printer layout information does not match M row dimension")
			for (j in 1:narrays) {
				spots <- 1:nspots
				for (gridr in 1:ngr)
				for (gridc in 1:ngc) {
					y <- object$M[spots,j]
					x <- object$A[spots,j]
					w <- weights[spots,j]
					object$M[spots,j] <- loessFit(y,x,w,span=span,iterations=iterations)$residuals
					spots <- spots + nspots
				}
			}
		},
		composite = {
			if(is.null(layout)) stop("Layout argument not specified")
			if(is.null(controlspots)) stop("controlspots argument not specified")
			ntips <- layout$ngrid.r * layout$ngrid.c
			nspots <- layout$nspot.r * layout$nspot.c
			for (j in 1:narrays) {
				y <- object$M[,j]
				x <- object$A[,j]
				w <- weights[,j]
				f <- is.finite(y) & is.finite(x)
				if(!is.null(w)) f <- f & is.finite(w)
				y[!f] <- NA
				fit <- loess(y~x,weights=w,span=span,subset=controlspots,na.action=na.exclude,degree=0,surface="direct",family="symmetric",trace.hat="approximate",iterations=iterations)
				alpha <- global <- y
				global[f] <- predict(fit,newdata=x[f])
				alpha[f] <- (rank(x[f])-1) / sum(f)
				spots <- 1:nspots
				for (tip in 1:ntips) {
					y <- object$M[spots,j]
					x <- object$A[spots,j]
					w <- weights[spots,j]
					local <- loessFit(y,x,w,span=span,iterations=iterations)$fitted
					object$M[spots,j] <- object$M[spots,j] - alpha[spots]*global[spots]-(1-alpha[spots])*local
					spots <- spots + nspots
				}
			}
		},
		control = {
#			if(is.null(layout)) stop("Layout argument not specified")
			if(is.null(controlspots)) stop("controlspots argument not specified")
#			ntips <- layout$ngrid.r * layout$ngrid.c
#			nspots <- layout$nspot.r * layout$nspot.c
			for (j in 1:narrays) {
				y <- object$M[,j]
				x <- object$A[,j]
				w <- weights[,j]
				f <- is.finite(y) & is.finite(x)
				if(!is.null(w)) f <- f & is.finite(w)
				y[!f] <- NA
				fit <- loess(y~x,weights=w,span=span,subset=controlspots,na.action=na.exclude,degree=1,surface="direct",family="symmetric",trace.hat="approximate",iterations=iterations)
				y[f] <- y[f]-predict(fit,newdata=x[f])
				object$M[,j] <- y
			}
		},
		robustspline = {
#			if(is.null(layout)) stop("Layout argument not specified")
			for (j in 1:narrays)
				object$M[,j] <- normalizeRobustSpline(object$M[,j],object$A[,j],layout,df=df,method=robust)
		}
	)
	object
}

normalizeRobustSpline <- function(M,A,layout=NULL,df=5,method="M") {
#	Robust spline normalization
#	Gordon Smyth
#	27 April 2003.  Last revised 14 January 2015.

	if(!requireNamespace("MASS",quietly=TRUE)) stop("MASS package required but is not installed (or can't be loaded)")
	if(!requireNamespace("splines",quietly=TRUE)) stop("splines package required but is not installed (or can't be loaded)")
	if(is.null(layout)) {
		ngrids <- 1
		nspots <- length(M)
	} else {
		ngrids <- round(layout$ngrid.r * layout$ngrid.c)
		nspots <- round(layout$nspot.r * layout$nspot.c)
		if(ngrids < 1) stop("layout incorrectly specified")
		if(nspots < df+1) stop("too few spots in each print-tip block")
	}

#	Global splines
	O <- is.finite(M) & is.finite(A)
	X <- matrix(NA,ngrids*nspots,df)
	X[O,] <- splines::ns(A[O],df=df,intercept=TRUE)
	x <- X[O,,drop=FALSE]
	y <- M[O]
	w <- rep(1,length(y))
	s0 <- summary(MASS::rlm(x,y,weights=w,method=method),method="XtWX",correlation=FALSE)
	beta0 <- s0$coefficients[,1]
	covbeta0 <- s0$cov * s0$stddev^2

#	Case with only one tip-group
	if(ngrids <= 1) {
		M[O] <- s0$residuals
		M[!O] <- NA
		return(M)
	}

#	Tip-wise splines
	beta <- array(1,c(ngrids,1)) %*% array(beta0,c(1,df))
	covbeta <- array(0,c(ngrids,df,df))
	spots <- 1:nspots
	for (i in 1:ngrids) {
		o <- O[spots]
		y <- M[spots][o]
		if(length(y)) {
			x <- X[spots,][o,,drop=FALSE]
			r <- qr(x)$rank
			if(r<df) x <- x[,1:r,drop=FALSE]
			w <- rep(1,length(y))
			s <- summary(MASS::rlm(x,y,weights=w,method=method),method="XtWX",correlation=FALSE)
			beta[i,1:r] <- s$coefficients[,1]
			covbeta[i,1:r,1:r] <- s$cov * s$stddev^2
		}
		spots <- spots + nspots
	}

#	Empirical Bayes estimates
	res.cov <- cov(beta) - apply(covbeta,c(2,3),mean)
	a <- max(0, sum(diag(res.cov)) / sum(diag(covbeta0)) )
	Sigma0 <- covbeta0 * a 
#	Sigma0 <- covbeta0 * max(0,mean(eigen(solve(covbeta0,res.cov))$values))

#	Shrunk splines
	if(a==0) {
		M[O] <- s0$residuals
		M[!O] <- NA
	} else {
		spots <- 1:nspots
		for (i in 1:ngrids) {
			beta[i,] <- beta0 + Sigma0 %*% solve(Sigma0+covbeta[i,,],beta[i,]-beta0)
			o <- O[spots]
			x <- X[spots,][o,]
			M[spots][o] <- M[spots][o] - x %*% beta[i,]
			M[spots][!o] <- NA
			spots <- spots + nspots
		}
	}
	M
}


#  PRINTORDER

normalizeForPrintorder <- function(object,layout,start="topleft",method="loess",separate.channels=FALSE,span=0.1,plate.size=32) {
#	Pre-normalize the foreground intensities for print order
#	Gordon Smyth
#	11 Mar 2002.  Last revised 18 June 2003.

	if(is.null(object$R) || is.null(object$G)) stop("List must contain components R and G")
	start <- match.arg(start,c("topleft","topright"))
	method <- match.arg(method,c("loess","plate"))
	po <- printorder(layout,start=start)$printorder
	nslides <- NCOL(object$R)
	for (i in 1:nslides) {
		RG <- normalizeForPrintorder.rg(R=object$R[,i],G=object$G[,i],printorder=po,method=method,separate.channels=separate.channels,span=span,plate.size=plate.size)
		object$R[,i] <- RG$R
		object$G[,i] <- RG$G
	}
	object
}

normalizeForPrintorder.rg <- function(R,G,printorder,method="loess",separate.channels=FALSE,span=0.1,plate.size=32,plot=FALSE) {
#	Pre-normalize the foreground intensities for print order, given R and G for a single array.
#	Gordon Smyth
#	8 Mar 2002.  Last revised 3 January 2007.

	if(plot) ord <- order(printorder)
	Rf <- log(R,2)
	Gf <- log(G,2)
	Rf[is.infinite(Rf)] <- NA
	Gf[is.infinite(Gf)] <- NA
	if(!separate.channels) Af <- (Rf+Gf)/2
	method <- match.arg(method,c("loess","plate"))
	if(method=="plate") {
		# Correct for plate pack (usually four 384-well plates)
		plate <- 1 + (printorder-0.5) %/% plate.size
		if(!requireNamespace("MASS",quietly=TRUE)) stop("MASS package required but is not installed (or can't be loaded)")
		hubermu <- function(...) MASS::huber(...)$mu
		if(separate.channels) {
			plate.mR <- tapply(Rf,plate,hubermu)
			plate.mG <- tapply(Gf,plate,hubermu)
			mR <- mG <- Rf
			for (i in 1:max(plate)) {
				mR[plate==i] <- plate.mR[i]
				mG[plate==i] <- plate.mG[i]
			}
			if(plot) {
				plot(printorder,Rf,xlab="Print Order",ylab="Log Intensity",type="n")
				points(printorder,Rf,pch=".",col="red")
				points(printorder,Gf,pch=".",col="green")
				lines(printorder[ord],mR[ord],col="red")
				lines(printorder[ord],mG[ord],col="green")
				return(invisible())
			}
			mR <- mR - mean(mR,na.rm=TRUE)
			mG <- mG - mean(mG,na.rm=TRUE)
		} else {
			plate.m <- tapply(Af,plate,hubermu)
			m <- Af
			for (i in 1:max(plate)) m[plate==i] <- plate.m[i]
			if(plot) {
				plot(printorder,Af,xlab="Print Order",ylab="Log Intensity",pch=".")
				lines(printorder[ord],m[ord])
			}
			mR <- mG <- m - mean(m,na.rm=TRUE)
		}
	} else {
		# Smooth correction for time order
		if(separate.channels) {
			mR <- fitted(loess(Rf~printorder,span=span,degree=0,family="symmetric",trace.hat="approximate",iterations=5,surface="direct",na.action=na.exclude))
			mG <- fitted(loess(Gf~printorder,span=span,degree=0,family="symmetric",trace.hat="approximate",iterations=5,surface="direct",na.action=na.exclude))
			if(plot) {
				plot(printorder,Rf,xlab="Print Order",ylab="Log Intensity",type="n")
				points(printorder,Rf,pch=".",col="red")
				points(printorder,Gf,pch=".",col="green")
				lines(printorder[ord],mR[ord],col="red")
				lines(printorder[ord],mG[ord],col="green")
				return(invisible())
			}
			mR <- mR - mean(mR,na.rm=TRUE)
			mG <- mG - mean(mG,na.rm=TRUE)
		} else {
			m <- fitted(loess(Af~printorder,span=span,degree=0,family="symmetric",trace.hat="approximate",iterations=5,surface="direct",na.action=na.exclude))
			if(plot) {
				plot(printorder,Af,xlab="Print Order",ylab="Log Intensity",pch=".")
				lines(printorder[ord],m[ord])
			}
			mR <- mG <- m - mean(m,na.rm=TRUE)
		}
	}
	list(R=2^(Rf-mR),G=2^(Gf-mG),R.trend=mR,G.trend=mG)
}

plotPrintorder <- function(object,layout,start="topleft",slide=1,method="loess",separate.channels=FALSE,span=0.1,plate.size=32) {
#	Pre-normalize the foreground intensities for print order.
#	Gordon Smyth
#	9 Apr 2002.  Last revised 18 June 2003.

	if(is.null(object$R) || is.null(object$G)) stop("List must contain components R and G")
	G <- object$G[,slide]
	R <- object$R[,slide]
	if(length(R) != length(G)) stop("R and G must have same length")
	start <- match.arg(start,c("topleft","topright"))
	po <- printorder(layout,start=start)$printorder
	invisible(normalizeForPrintorder.rg(R=R,G=G,printorder=po,method=method,separate.channels=separate.channels,span=span,plate.size=plate.size,plot=TRUE))
}

#  BETWEEN ARRAY NORMALIZATION

normalizeBetweenArrays <- function(object, method=NULL, targets=NULL, cyclic.method="fast", ...) {
#	Normalize between arrays
#	Gordon Smyth
#	12 Apri 2003.  Last revised 13 June 2016.

	if(is.data.frame(object)) {
		object <- as.matrix(object)
		if(mode(object) != "numeric") stop("'object' is a data.frame and not all columns are numeric")
	}

#	Default method
	if(is.null(method)) {
		if(is(object,"matrix")) {
			method="quantile"
		} else if(is(object,"EListRaw")) {
			method="quantile"
		} else {
			method="Aquantile"
		}	
	}
	choices <- c("none","scale","quantile","Aquantile","Gquantile","Rquantile","Tquantile","vsn","cyclicloess")
	method <- match.arg(method,choices)
	if(method=="vsn") stop("vsn method no longer supported. Please use normalizeVSN instead.")

#	Methods for matrices
	if(is(object,"matrix")) {
		if(!(method %in% c("none","scale","quantile","cyclicloess"))) stop("method not applicable to single-channel data")
		return(switch(method,
			none = object,
			scale = normalizeMedianValues(object),
			quantile = normalizeQuantiles(object, ...),
			cyclicloess = normalizeCyclicLoess(object,method=cyclic.method, ...)
		))
	}

#	Treat EListRaw objects as matrices
	if(is(object,"EListRaw")) {
		if(method=="cyclicloess")
			object$E <- Recall(log2(object$E),method=method,...)
		else
			object$E <- log2(Recall(object$E,method=method,...))
		object <- new("EList",unclass(object))
		return(object)
	}

#	From here assume two-color data
	if(is(object,"RGList")) object <- MA.RG(object)
	if(is.null(object$M) || is.null(object$A)) stop("object doesn't appear to be RGList or MAList object")
	switch(method,
		scale = {
			object$M <- normalizeMedianAbsValues(object$M)
			object$A <- normalizeMedianAbsValues(object$A)
		},
		quantile = {
			narrays <- NCOL(object$M)
			Z <- normalizeQuantiles(cbind(object$A-object$M/2,object$A+object$M/2),...)
			G <- Z[,1:narrays]
			R <- Z[,narrays+(1:narrays)]
			object$M <- R-G
			object$A <- (R+G)/2
		},
		Aquantile = {
			object$A <- normalizeQuantiles(object$A,...)
		},
		Gquantile = {
			G <- object$A-object$M/2
			E <- normalizeQuantiles(G,...) - G
			object$A <- object$A + E
		},
		Rquantile = {
			R <- object$A+object$M/2
			E <- normalizeQuantiles(R,...) - R
			object$A <- object$A + E
		},
		Tquantile = {
			narrays <- NCOL(object$M)
			if(NCOL(targets)>2) targets <- targets[,c("Cy3","Cy5")]
			targets <- as.vector(targets)
			Z <- cbind(object$A-object$M/2,object$A+object$M/2)
			for (u in unique(targets)) {
				j <- targets==u
				Z[,j] <- normalizeQuantiles(Z[,j],...)
			}
			G <- Z[,1:narrays]
			R <- Z[,narrays+(1:narrays)]
			object$M <- R-G
			object$A <- (R+G)/2
		})
	object
}


normalizeVSN <- function(x,...)
{
	if(!requireNamespace("Biobase",quietly=TRUE)) stop("Biobase package required but is not installed (or can't be loaded)")
	if(!requireNamespace("vsn",quietly=TRUE)) stop("vsn package required but is not installed (or can't be loaded)")
	UseMethod("normalizeVSN")
}

normalizeVSN.RGList <- function(x,...)
#	vsn background correction and normalization for RGList objects
#	Gordon Smyth
#	9 Sep 2010.
{
	x <- backgroundCorrect(x,method="subtract")
	y <- cbind(x$G,x$R)
	x$G <- x$R <- NULL
	y <- Biobase::exprs(vsn::vsnMatrix(x=y,...))
	n2 <- ncol(y)/2
	G <- y[,1:n2]
	R <- y[,n2+(1:n2)]
	x$M <- R-G
	x$A <- (R+G)/2
	new("MAList",unclass(x))
}

normalizeVSN.EListRaw <- function(x,...)
#	vsn background correction and normalization for EListRaw objects
#	Gordon Smyth
#	9 Sep 2010.
{
	x <- backgroundCorrect(x,method="subtract")
	x$E <- Biobase::exprs(vsn::vsnMatrix(x=x$E,...))
	new("EList",unclass(x))
}

normalizeVSN.default <- function(x,...)
#	vsn background correction and normalization for matrices
#	Gordon Smyth
#	9 Sep 2010.
{
	Biobase::exprs(vsn::vsnMatrix(x=as.matrix(x),...))
}


normalizeQuantiles <- function(A, ties=TRUE) {
#	Normalize columns of a matrix to have the same quantiles, allowing for missing values.
#	Gordon Smyth
#	25 June 2002.  Last revised 16 May 2019.

	n <- dim(A)
	if(is.null(n)) return(A)
	if(n[2]==1) return(A)
	O <- S <- array(,n)
	nobs <- rep(n[1],n[2])
	i <- (0:(n[1]-1))/(n[1]-1)
	for (j in 1:n[2]) {
		Si <- sort(A[,j], method="quick", index.return=TRUE)
		nobsj <- length(Si$x)
		if(nobsj < n[1]) {
			nobs[j] <- nobsj
			isna <- is.na(A[,j])
			S[,j] <- approx((0:(nobsj-1))/(nobsj-1), Si$x, i, ties=list("ordered",mean))$y
			O[!isna,j] <- ((1:n[1])[!isna])[Si$ix]
		} else {
			S[,j] <- Si$x
			O[,j] <- Si$ix
		}
	}
	m <- rowMeans(S)
	for (j in 1:n[2]) {
		if(ties) r <- rank(A[,j])
		if(nobs[j] < n[1]) {
			isna <- is.na(A[,j])
			if(ties)
				A[!isna,j] <- approx(i, m, (r[!isna]-1)/(nobs[j]-1), ties=list("ordered",mean))$y
			else
				A[O[!isna,j],j] <- approx(i, m, (0:(nobs[j]-1))/(nobs[j]-1), ties=list("ordered",mean))$y
		} else {
			if(ties)
				A[,j] <- approx(i, m, (r-1)/(n[1]-1), ties=list("ordered",mean))$y
			else
				A[O[,j],j] <- m
		}
	}
	A
}

normalizeMedianAbsValues <- function(x) 
#	Normalize columns of a matrix to have the same median absolute value
#	Gordon Smyth
#	12 April 2003.  Last modified 19 Oct 2005.
{
	narrays <- NCOL(x)
	if(narrays==1) return(x)
	cmed <- log(apply(abs(x), 2, median, na.rm=TRUE))
	cmed <- exp(cmed - mean(cmed))
	t(t(x)/cmed)
}

normalizeMedianValues <- function(x) 
#	Normalize columns of a matrix to have the same median value
#	Gordon Smyth
#	24 Jan 2011.  Last modified 24 Jan 2011.
{
	narrays <- NCOL(x)
	if(narrays==1) return(x)
	cmed <- log(apply(x, 2, median, na.rm=TRUE))
	cmed <- exp(cmed - mean(cmed))
	t(t(x)/cmed)
}

normalizeCyclicLoess <- function(x, weights = NULL, span=0.7, iterations = 3, method="fast")
#	Cyclic loess normalization of columns of matrix
#	incorporating probes weights.
#	Yunshun (Andy) Chen and Gordon Smyth
#	14 April 2010.  Last modified 24 Feb 2012.
{
	x <- as.matrix(x)
	method <- match.arg(method, c("fast","affy","pairs"))
	n <- ncol(x)
	if(method=="pairs") {
		for (k in 1:iterations)
		for (i in 1:(n-1)) for (j in (i+1):n) {
			m <- x[,j] - x[,i]
			a <- .5*(x[,j] + x[,i])
			f <- loessFit(m, a, weights=weights, span=span)$fitted
			x[,i] <- x[,i] + f/2
			x[,j] <- x[,j] - f/2		
		}
	}
	if(method=="fast") {
		for (k in 1:iterations) {
			a <- rowMeans(x,na.rm=TRUE)
			for (i in 1:n){
				m <- x[,i] - a
				f <- loessFit(m, a, weights=weights, span=span)$fitted
				x[,i] <- x[,i] - f
			}
		}
	}
	if(method=="affy") {
		g <- nrow(x)
		for (k in 1:iterations) {
			adjustment <- matrix(0,g,n)
			for (i in 1:(n-1)) for (j in (i+1):n) {
				m <- x[,j] - x[,i]
				a <- .5*(x[,j] + x[,i])
				f <- loessFit(m, a, weights=weights, span=span)$fitted
				adjustment[,j] <- adjustment[,j] + f
				adjustment[,i] <- adjustment[,i] - f
			}
			x <- x-adjustment/n
		}
	}
	x
}


