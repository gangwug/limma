fitmixture <- function(log2e,mixprop,niter=4,trace=FALSE)
#	Fit mixture model by non-linear least squares
#	Gordon Smyth  9 April 2007
{
	log2e <- as.matrix(log2e)
	mixprop <- as.numeric(mixprop)
	narrays <- ncol(log2e)
	nprobes <- nrow(log2e)
	y <- 2^log2e
	start <- lm.fit(cbind(mixprop,1-mixprop),t(y))$coef
	start <- pmax(start,1)
	b <- as.vector(c(1,-1) %*% log(start))
	z <- log(y)
	pm <- matrix(1,nprobes,1) %*% matrix(2*mixprop-1,1,narrays)
	for (i in 1:niter) {
		mub <- logcosh(b/2)+log(1+tanh(b/2)*pm)
		a <- rowMeans(z-mub)
		mu <- a+mub
		if(trace) {
			s <- sqrt(narrays/(narrays-2)*rowMeans((z-mu)^2))
			if(i>1) cat("stdev changes",summary(sold-s),"\n")
			sold <- s
		}
		tb <- tanh(b/2)
		dmu <- (tb+pm)/(1+tb*pm)/2
		b <- b + rowMeans(dmu*(z-mu)) / (1e-6+rowMeans((dmu-rowMeans(dmu))^2))
		if(trace) cat("betas",summary(b),"\n")
	}
	mub <- logcosh(b/2)+log(1+tanh(b/2)*pm)
	a <- rowMeans(z-mub)
	mu <- a+mub
	s <- sqrt(narrays/(narrays-2)*rowMeans((z-mu)^2))
	l2 <- log(2)
	list(A=a/l2,M=b/l2,stdev=s/l2)
}

