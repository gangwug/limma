fitGammaIntercept <- function(y,offset=0,maxit=1000)
#	Estimate intercept of additive gamma glm
#	Belinda Phipson and Gordon Smyth
#	Created 17 May 2012.  Last modified 2 Dec 2012.
{
	if(any(y<0)) stop("negative y not permitted")
	if(any(offset<0)) stop("offsets must be positive")

#	Treat constant offset as special case
	r <- range(offset)
	if(r[1]+1e-14 > r[2]) return(mean(y)-r[1])

	n <- length(y)
	x <- 0
	iter <- 0
	repeat {
		iter <- iter+1
		Q <- sum(y/(offset+x))
		dQ <- sum(y/(offset+x)^2)
		dif <- (Q-n)/dQ
		x <- x+dif
		if(abs(dif) < 1e-8) break
		if(iter > maxit) {
			warning("iteration limit maxit exceeded")
			break
		}
	}
	x
}
