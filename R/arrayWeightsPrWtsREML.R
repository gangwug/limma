.arrayWeightsPrWtsREML <- function(y, design=NULL, weights, var.design, prior.n=10, maxiter=50L, tol=1e-6, trace=FALSE)
#	Estimate array weights by REML allowing for prior observation weights
#	Probes with missing or infinite values are removed.
#	Uses an exact Fisher scoring algorithm similar to statmod::remlscor.
#	Gordon Smyth
#	Created 12 Feb 2019 from .arrayWeightsREML.
#	Last revised 15 Feb 2019.
{
#	y should be a numeric matrix
	narrays <- ncol(y)
	ngenes <- nrow(y)
	p <- ncol(design)

#	Columns of var.design should sum to zero, and intercept column should be omitted.
	Z2 <- var.design
	Z <- cbind(1,Z2)
	ngam <- ncol(Z2)

#	Starting values
	iter <- 0L
	gam <- rep_len(0,ngam)
	w <- rep_len(1,narrays)
	if(trace) cat("iter convcrit range(w)\n")

#	Fisher scoring iteration
	p2 <- (p * (p+1L)) %/% 2L
	p1p2 <- (p+1L):p2
	Eye <- diag(1,nrow=narrays,ncol=p)
	sqrt2 <- sqrt(2)
	info2 <- matrix(0,ngam,ngam)
	repeat {
		iter <- iter+1L

#		For shrinkage of gam, start from prior weight of prior.n genes
		info2 <- prior.n*crossprod(Z2)
		z <- prior.n*(w-1)

		for (g in 1:ngenes) {
#			Fit weighted linear model and extract residual variances
			fitm <- lm.wfit(design, y[g,], w*weights[g,])
			s2 <- mean(fitm$effects[(fitm$rank+1):narrays]^2)

#			Fisher information matrix for variance parameters (including intercept)
			Q <- qr.qy(fitm$qr,Eye)
			Q2 <- matrix(0,narrays,p2)
			j0 <- 0L
			for (k in 0:(p-1L)) {
				Q2[,(j0+1L):(j0+p-k)] <- Q[,1:(p-k)]*Q[,(k+1):p]
				j0 <- j0+p-k
			}
			if(p > 1) Q2[,p1p2] <- sqrt2*Q2[,p1p2]
			h <- rowSums(Q2[,1:p,drop=FALSE])
			info <- crossprod(Z,(1-2*h)*Z) + crossprod(crossprod(Q2,Z))

#			Fisher information excluding intercept (i.e., for gam)
			info2 <- info2 + info[-1,-1,drop=FALSE] - (info[-1,1,drop=FALSE]/info[1,1]) %*% info[1,-1,drop=FALSE]

#			Variance model residual
			if(s2 > 1e-15) {
				z <- z + w * weights[g,] * fitm$residuals^2 / s2 - (1-h)
			}
		}

#		Average information matrix and variance residual per gene
		info2 <- info2 / (ngenes+prior.n)
		z <- z / (ngenes+prior.n)

#		Fisher scoring
		dl <- crossprod(Z2, z)
		gamstep <- solve(info2, dl)
		gam <- gam + gamstep

#		Update array weights
		w <- drop(exp(Z2 %*% (-gam)))

#		Test for convergence
		convcrit <- crossprod(dl,gamstep) / (ngenes+prior.n) / ngam
		if(trace) cat(iter,convcrit,range(w),"\n")
		if(is.na(convcrit)) {
			warning("convergence tolerance not achievable, stopping prematurely")
			break
		}
		if(convcrit < tol) break

		if(iter==maxiter) {
			warning("iteration limit reached")
			break
		}
	}
	w
}
