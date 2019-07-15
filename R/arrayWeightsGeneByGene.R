.arrayWeightsGeneByGene <- function(E, design=NULL, weights=NULL, var.design=NULL, prior.n=10, trace=FALSE)
#	Estimate array variances via gene-by-gene update algorithm
#	Created by Matt Ritchie 7 Feb 2005.
#	Cynthia Liu added var.design argument 22 Sep 2014.
#	Fixes and speedups by Gordon Smyth 15 Feb 2019.
#	Last modified 14 Feb 2019.
{
	ngenes <- nrow(E)
	narrays <- ncol(E)
	if(is.null(design)) design <- matrix(1,narrays,1)
	nparams <- ncol(design)

#	Columns of var.design should sum to zero
	if(is.null(var.design)) {
		Z2 <- contr.sum(narrays)
	} else {
		Z2 <- var.design
	}
	Z <- cbind(1,Z2)

#	Intialise array gammas to zero (with prior weight of prior.n genes having leverage=0)
	ngam <- ncol(Z2)
	gam <- rep_len(0, ngam)
	aw <- rep_len(1,narrays)
	info2 <- prior.n*crossprod(Z2)

#	If requested, progess will be output 10 times at equal intervals
	if(trace) {
		cat("gene range(w)\n")
		ReportInterval <- pmax(as.integer(ngenes/10),1L)
	}

#	Step progressive algorithm once through all genes
	Zero <- rep_len(0,narrays)
	One <- rep_len(1,narrays)
	for(i in 1:ngenes) {
		if(is.null(weights)) {
			w <- aw
		} else {
			w <- aw*weights[i,]
		}
		y <- E[i,]
		if(anyNA(y)) {
			obs <- is.finite(y)
			nobs <- sum(obs)
			if(nobs <= 2L) next
			X <- design[obs, , drop = FALSE]
			y <- y[obs]
			w <- w[obs]
			fit <- lm.wfit(X, y, w)
			h1 <- d <- Zero
			d[obs] <- w*fit$residuals^2
			h1[obs] <- 1-hat(fit$qr)
		} else {
			fit <- lm.wfit(design, y, w)
			d <- w*fit$residuals^2
			h1 <- 1-hat(fit$qr)
		}
		s2 <- mean(fit$effects[-(1:fit$rank)]^2)
		if(s2 < 1e-15) next
		info <- crossprod(Z, h1*Z)
		info2 <- info2 + info[-1,-1,drop=FALSE] - (info[-1,1,drop=FALSE]/info[1,1]) %*% info[1,-1,drop=FALSE]
		z <- d/s2 - h1
		dl <- crossprod(Z2, z)
		gam <- gam + solve(info2, dl)
	 	aw <- drop(exp(Z2 %*% (-gam)))

#		Progress output
		if(trace && (i %% ReportInterval==0L)) cat(i,range(aw),"\n")
	}

	aw
}
