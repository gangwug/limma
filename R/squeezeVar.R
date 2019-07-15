#	EMPIRICAL BAYES SQUEEZING OF VARIANCES

squeezeVar <- function(var, df, covariate=NULL, robust=FALSE, winsor.tail.p=c(0.05,0.1))
#	Empirical Bayes posterior variances
#	Gordon Smyth
#	Created 2 March 2004.  Last modified 5 May 2016.
{
	n <- length(var)

#	Degenerate special cases
	if(n == 0) stop("var is empty")
	if(n == 1) return(list(var.post=var,var.prior=var,df.prior=0))

#	When df==0, guard against missing or infinite values in var
	if(length(df)>1) var[df==0] <- 0

#	Estimate hyperparameters
	if(robust) {
		fit <- fitFDistRobustly(var, df1=df, covariate=covariate, winsor.tail.p=winsor.tail.p)
		df.prior <- fit$df2.shrunk
	} else {
		fit <- fitFDist(var, df1=df, covariate=covariate)
		df.prior <- fit$df2
	}
	if(anyNA(df.prior)) stop("Could not estimate prior df")

#	Posterior variances
	var.post <- .squeezeVar(var=var, df=df, var.prior=fit$scale, df.prior=df.prior)

	list(df.prior=df.prior,var.prior=fit$scale,var.post=var.post)
}

.squeezeVar <- function(var, df, var.prior, df.prior)
#	Squeeze posterior variances given hyperparameters
#	NAs not allowed in df.prior
#	Gordon Smyth
#	Created 5 May 2016
{
	n <- length(var)
	isfin <- is.finite(df.prior)
	if(all(isfin)) return( (df*var + df.prior*var.prior) / (df+df.prior) )

#	From here, at least some df.prior are infinite

#	For infinite df.prior, return var.prior
	if(length(var.prior) == n) {
		var.post <- var.prior
	} else {
		var.post <- rep_len(var.prior, length.out=n)
	}

#	Maybe some df.prior are finite
	if(any(isfin)) {
		i <- which(isfin)
		if(length(df)>1) df <- df[i]
		df.prior <- df.prior[i]
		var.post[i] <- (df*var[i] + df.prior*var.post[i]) / (df+df.prior)
	}

	var.post
}
