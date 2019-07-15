predFCm <- function(fit,coef=2, var.indep.of.fc=TRUE, all.de=TRUE, prop.true.null.method="lfdr")
#	Predictive (empirical Bayes shrunk) fold changes
#	Belinda Phipson and Gordon Smyth
#	Created 29 May 2012. Last modified 1 December 2014.
{
#	Check fit
	if(!is(fit,"MArrayLM")) stop("fit must be a MArrayLM fit object")
	if(is.null(fit$p.value)) fit <- eBayes(fit)

#	Check coef
	coef <- coef[1]

#	Estimate proportion of true nulls and re-run eBayes
	p <- 1-propTrueNull(fit$p.value[,coef], method=prop.true.null.method)
	if(p==0) p <- 1e-8
	if(length(fit$s2.prior)==1L) trend <- FALSE else trend <- TRUE
	if(length(fit$df.prior)==1L) robust <- FALSE else robust <- TRUE
	fit <- eBayes(fit,proportion=p,trend=trend,robust=robust)
	v <- fit$cov.coefficients[coef,coef]

	if(var.indep.of.fc){
		v0 <- fitGammaIntercept(fit$coeff[,coef]^2, offset=v*fit$s2.post)
		if(v0<0) v0<-1e-8
		pfc <- fit$coeff[,coef]*v0/(v0+v*fit$s2.post)
		if(!all.de) {
			A <- p/(1-p)
			B <- (v*fit$s2.post/(v*fit$s2.post+v0))^0.5
			C <- exp(fit$coeff[,coef]^2*v0/(2*v^2*fit$s2.post^2+2*v*v0*fit$s2.post))
			lods <- log(A*B*C)
			probDE <- exp(lods)/(1+exp(lods))
			probDE[lods>700] <- 1
			pfc <- pfc*probDE
		}
	} else {
		if(is.null(fit$s2.post)) fit <- eBayes(fit)
		b2 <- fit$coeff[,coef]^2 / fit$s2.post
		v0 <- fitGammaIntercept(b2, offset=v)
		v0 <- pmin(v0, 1e-8)
		pfc <- fit$coeff[,coef] * v0/(v0+v)
		if(!all.de) {
			A <- p/(1-p)
			B <- (v/(v+v0))^0.5
			C <- exp(fit$coeff[,coef]^2*v0/(2*v^2*fit$s2.post+2*v*v0*fit$s2.post))
			lods <- log(A*B*C)
			probDE <- exp(lods)/(1+exp(lods))
			probDE[lods>700] <- 1
			pfc <- pfc*probDE
		}
	}
	pfc
}
