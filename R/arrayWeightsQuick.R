arrayWeightsQuick <- function(y, fit)
#	Compute approximate array quality weights from a MArrayLM object
#	Gordon Smyth
#	25 Oct 2004.  Last revised 28 Oct 2004.
{
	if(!is.null(fit$weights)) warning("spot quality weights found but not taken into account")
	res <- as.matrix(y)- fit$coef %*% t(fit$design)
	h <- hat(fit$design, intercept=FALSE)
	mures2 <- fit$sigma^2 %*% array(1-h,c(1,length(h)))
	1/colMeans(res*res/mures2,na.rm=TRUE)
}
