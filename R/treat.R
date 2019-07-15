###  treat.R

treat <- function(fit, lfc=log2(1.2), trend=FALSE, robust=FALSE, winsor.tail.p=c(0.05,0.1))
#	Moderated t-statistics with threshold
#	Davis McCarthy, Gordon Smyth
#	25 July 2008.  Last revised 25 March 2018.
{
#	Check fit
	if(!is(fit,"MArrayLM")) stop("fit must be an MArrayLM object")
	if(is.null(fit$coefficients)) stop("coefficients not found in fit object")
	if(is.null(fit$stdev.unscaled)) stop("stdev.unscaled not found in fit object")
	fit$lods <- NULL

	coefficients <- as.matrix(fit$coefficients)
	stdev.unscaled <- as.matrix(fit$stdev.unscaled)
	sigma <- fit$sigma
	df.residual <- fit$df.residual
	if (is.null(coefficients) || is.null(stdev.unscaled) || is.null(sigma) || 
		is.null(df.residual)) 
		stop("No data, or argument is not a valid lmFit object")
	if (all(df.residual == 0)) 
		stop("No residual degrees of freedom in linear model fits")
	if (all(!is.finite(sigma))) 
		stop("No finite residual standard deviations")
	if(trend) {
		covariate <- fit$Amean
		if(is.null(covariate)) stop("Need Amean component in fit to estimate trend")
	} else {
		covariate <- NULL
	}
	sv <- squeezeVar(sigma^2, df.residual, covariate=covariate, robust=robust, winsor.tail.p=winsor.tail.p)
	fit$df.prior <- sv$df.prior
	fit$s2.prior <- sv$var.prior
	fit$s2.post <- sv$var.post
	df.total <- df.residual + sv$df.prior
	df.pooled <- sum(df.residual,na.rm=TRUE)
	df.total <- pmin(df.total,df.pooled)
	fit$df.total <- df.total
	lfc <- abs(lfc)
	acoef <- abs(coefficients)
	se <- stdev.unscaled*sqrt(fit$s2.post)
	tstat.right <- (acoef-lfc)/se
	tstat.left <- (acoef+lfc)/se
	fit$t <- array(0,dim(coefficients),dimnames=dimnames(coefficients))
	fit$p.value <- pt(tstat.right, df=df.total,lower.tail=FALSE) + pt(tstat.left,df=df.total,lower.tail=FALSE)
	tstat.right <- pmax(tstat.right,0)
	fc.up <- (coefficients > lfc)
	fc.down <- (coefficients < -lfc)
	fit$t[fc.up] <- tstat.right[fc.up]
	fit$t[fc.down] <- -tstat.right[fc.down]
	fit$treat.lfc <- lfc
	fit
}

topTreat <- function(fit,coef=1,sort.by="p",resort.by=NULL,...)
#	Summary table of top genes by treat
#	Gordon Smyth
#	15 June 2009.  Last modified 27 February 2014.

#	NOTE: This function may become deprecated soon as topTable() takes over
{
#	Check coef is length 1
	if(length(coef)>1) {
		coef <- coef[1]
		warning("Treat is for single coefficients: only first value of coef being used")
	}

#	Check sort.by and resort.by
	if(sort.by=="B") stop("Trying to sort.by B, but treat doesn't produce a B-statistic")
	if(!is.null(resort.by)) if(resort.by=="B") stop("Trying to resort.by B, but treat doesn't produce a B-statistic")

	topTable(fit=fit,coef=coef,sort.by=sort.by,...)
}

