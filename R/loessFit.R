#	LOESS FUNCTIONS

loessFit <- function(y, x, weights=NULL, span=0.3, iterations=4L, min.weight=1e-5, max.weight=1e5, equal.weights.as.null=TRUE, method="weightedLowess")
#	Fast lowess fit for univariate x and y allowing for weights
#	Uses lowess() if weights=NULL and weightedLowess(), locfit.raw() or loess() otherwise
#	Gordon Smyth
#	28 June 2003.  Last revised 14 January 2015.
{
#	Check x and y
	n <- length(y)
	if(length(x) != n) stop("y and x have different lengths")
	out <- list(fitted=rep(NA,n),residuals=rep(NA,n))

	obs <- is.finite(y) & is.finite(x)
	xobs <- x[obs]
	yobs <- y[obs]
	nobs <- length(yobs)

#	If no good obs, exit straight away
	if(nobs==0) return(out)

#	Check span
	if(span < 1/nobs) {
		out$fitted[obs] <- y[obs]
		out$residuals[obs] <- 0
		return(out)
	}

#	Check min.weight
	if(min.weight<0) min.weight <- 0

#	Check weights
	if(!is.null(weights)) {
		if(length(weights) != n) stop("y and weights have different lengths")
		wobs <- weights[obs]
		wobs[is.na(wobs)] <- 0
		wobs <- pmax(wobs,min.weight)
		wobs <- pmin(wobs,max.weight)
#		If weights all equal, treat as NULL
		if(equal.weights.as.null) {
			r <- range(wobs)
			if(r[2]-r[1] < 1e-15) weights <- NULL
		}
	}

#	If no weights, so use classic lowess algorithm
	if(is.null(weights)) {
		o <- order(xobs)
		lo <- lowess(x=xobs,y=yobs,f=span,iter=iterations-1L)
		out$fitted[obs][o] <- lo$y
		out$residuals[obs] <- yobs-out$fitted[obs]
		return(out)
	}

#	Count number of observations with positive weights (must always be positive)
	if(min.weight>0)
		nwobs <- nobs
	else 
		nwobs <- sum(wobs>0)

#	Check whether too few obs to estimate lowess curve
	if(nwobs < 4+1/span) {
		if(nwobs==1L) {
			out$fitted[obs] <- yobs[wobs>0]
			out$residuals[obs] <- yobs-out$fitted[obs]
		} else {
			fit <- lm.wfit(cbind(1,xobs),yobs,wobs)
			out$fitted[obs] <- fit$fitted
			out$residuals[obs] <- fit$residuals
		}
		return(out)
	}

#	Need to compute lowess with unequal weights
	method <- match.arg(method, c("weightedLowess","locfit","loess"))
	switch(method,
		"weightedLowess" = {
			fit <- weightedLowess(x=xobs,y=yobs,weights=wobs,span=span,iterations=iterations,npts=200)
			out$fitted[obs] <- fit$fitted
			out$residuals[obs] <- fit$residuals
		},
		"locfit" = {
#			Check for locfit package
			if(!requireNamespace("locfit",quietly=TRUE)) stop("locfit required but is not installed (or can't be loaded)")
#			Weighted lowess with robustifying iterations
		    biweights <- rep(1,nobs)
 			for (i in 1:iterations) {
       			fit <- locfit::locfit(yobs~xobs,weights=wobs*biweights,alpha=span,deg=1)
		        res <- residuals(fit,type="raw")
		        s <- median(abs(res))
		        biweights <- pmax(1-(res/(6*s))^2,0)^2
		    }
			out$fitted[obs] <- fitted(fit)
			out$residuals[obs] <- res
		},
		"loess" = {
#			Suppress warning "k-d tree limited by memory"
			oldopt <- options(warn=-1)
			on.exit(options(oldopt))
			bin <- 0.01
			fit <- loess(yobs~xobs,weights=wobs,span=span,degree=1,parametric=FALSE,normalize=FALSE,statistics="approximate",surface="interpolate",cell=bin/span,iterations=iterations,trace.hat="approximate")
			out$fitted[obs] <- fit$fitted
			out$residuals[obs] <- fit$residuals
		}
	)
	out
}

