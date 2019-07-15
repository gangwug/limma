.arrayWeightsREML <- function(y,design=NULL,var.design=NULL,prior.n=10,maxiter=50L,tol=1e-6,trace=FALSE)
#	Estimate array weights by REML.
#	Assumes no observation weights, missing values or infinite values.
#	Uses an exact Fisher scoring algorithm similar to statmod::remlscor.
#	Gordon Smyth
#	Created 13 Dec 2005. Last revised 15 Feb 2019.
{
#	y should be a numeric matrix
	narrays <- ncol(y)
	ngenes <- nrow(y)
	if(is.null(design)) design <- matrix(1,narrays,1)
	p <- ncol(design)

#	Columns of var.design should sum to zero, and intercept column should be omitted.
	if(is.null(var.design)) {
		Z2 <- contr.sum(narrays)
	} else {
		Z2 <- var.design
	}
	Z <- cbind(1,Z2)
	ngam <- ncol(Z2)

#	Starting values
	iter <- 0L
	gam <- rep_len(0,ngam)
	w <- rep_len(1,narrays)
	convcrit.last <- Inf
	if(trace) cat("iter convcrit range(w)\n")

#	Initial fit of unweighted linear models to check for zero variances
	fitm <- lm.fit(design, t(y))
	Effects <- fitm$effects[(fitm$rank+1):narrays,,drop=FALSE]
	s2 <- colMeans(Effects^2)

#	Remove all rows with no residual variance
	if(min(s2) < 1e-15) {
		ok <- which(s2 >= 1e-15)
		y <- y[ok, ,drop=FALSE]
		ngenes <- nrow(y)
		if(ngenes < 2L) {
			names(w) <- colnames(y)
			return(w)
		}
		Effects <- Effects[,ok,drop=FALSE]
		fitm$residuals <- fitm$residuals[,ok,drop=FALSE]
		s2 <- s2[ok]
	}

#	Fisher scoring iteration
	p2 <- (p * (p+1L)) %/% 2L
	Q2 <- array(0,c(narrays,p2))
	repeat {
		iter <- iter+1L

#		Fit weighted linear models and extract residual variances
		if(iter > 1L) {
			fitm <- lm.wfit(design, t(y), w)
			Effects <- fitm$effects[(fitm$rank+1):narrays,,drop=FALSE]
			s2 <- colMeans(Effects^2)
		}

#		Fisher information matrix for variance parameters (including intercept)
		Q <- qr.qy(fitm$qr,diag(1,nrow=narrays,ncol=p))
		j0 <- 0
		for (k in 0:(p-1)) {
			Q2[,(j0+1):(j0+p-k)] <- Q[,1:(p-k)]*Q[,(k + 1):p]
			j0 <- j0 + p - k
		}
		if(p > 1) Q2[,(p+1):p2] <- sqrt(2)*Q2[,(p+1):p2]
		h <- rowSums(Q2[,1:p,drop=FALSE])
		info <- crossprod(Z,(1-2*h)*Z) + crossprod(crossprod(Q2,Z))

#		Fisher information excluding intercept (i.e., for gam)
		info2 <- info[-1,-1,drop=FALSE] - (info[-1,1,drop=FALSE]/info[1,1]) %*% info[1,-1,drop=FALSE]

#		Score vector (log-lik derivative) for gam
		z <- t(w * fitm$residuals^2) / s2
		z <- colMeans(z) - (1-h)

#		Add prior support for w=1 and gam=0
		info2 <- ngenes*info2 + prior.n*crossprod(Z2)
		z <- ngenes*z + prior.n*(w-1)

#		Fisher scoring
		dl <- crossprod(Z2, z)
		gamstep <- solve(info2,dl)

#		Convergence criterion
		convcrit <- crossprod(dl,gamstep) / ngam / (ngenes+prior.n)
		if(is.na(convcrit) || convcrit >= convcrit.last) {
			warning("convergence tolerance not achievable, stopping prematurely")
			break
		} else {
			convcrit.last <- convcrit
		}

#		Update array weights
		gam <- gam + gamstep
		w <- drop(exp(Z2 %*% (-gam)))
		if(trace) cat(iter,convcrit,range(w),"\n")
#		if(trace) cat("gamstep sd ",sqrt(sum(gamstep^2)),"\n")

#		Test for convergence
		if(convcrit < tol) break

#		Check for iter max
		if(iter==maxiter) {
			warning("iteration limit reached")
			break
		}
	}
	w
}
