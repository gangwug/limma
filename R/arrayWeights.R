arrayWeights <- function(object, design=NULL, weights=NULL, var.design=NULL, var.group=NULL, prior.n=10, method="auto", maxiter=50L, tol=1e-5, trace=FALSE)
#	Estimate array quality weights.
#
#	Created by Matt Ritchie 7 Feb 2005.
#	Gordon Smyth simplified argument checking to use getEAWP, 9 Mar 2008.
#	Cynthia Liu added var.design argument 22 Sep 2014.
#	Rewrite by Gordon Smyth 12 Feb 2019.
#	Last modified 14 Feb 2019.
{
#	Check object
	y <- getEAWP(object)
	E <- y$exprs
	ngenes <- nrow(E)
	narrays <- ncol(E)

#	Initial values for array weights
	w <- rep_len(1,narrays)
	names(w) <- colnames(E)

#	Require at least 2 rows for estimates to be useful
	if(ngenes < 2L) return(w)

#	Check design
	if(is.null(design)) {
		design <- matrix(1,narrays,1)
		p <- 1L
	} else {
		design <- as.matrix(design)
		if(mode(design) != "numeric") stop("design must be a numeric matrix")
		QR <- qr(design)
		p <- QR$rank

#		If not full rank, remove superfluous columns
		if(p < ncol(design)) design <- design[,QR$pivot[1:p],drop=FALSE]
	}

#	Require at least 2 residual df.
	if(narrays - p < 2L) return(w)

#	Check weights
	if(is.null(weights)) weights <- y$weights
	if(!is.null(weights)) {
		weights <- as.matrix(weights)
		if(nrow(weights)!=ngenes || ncol(weights)!=narrays) stop("dimensions of weights should match those of expression object")
		r <- range(weights)
		if(is.na(r[1])) stop("NA weights not allowed")
		if(r[1] < 0) stop("Negative weights not allowed")
		if(is.infinite(r[2])) stop("Infinite weights not allowed")
		if(r[1]==0) {
			E[weights==0] <- NA
			weights[weights==0] <- 1
		}
	}

#	var.group takes precedence over var.design
	if(!is.null(var.group)) {
		var.group <- droplevels(as.factor(var.group))
		if(length(var.group) != narrays) stop("var.group has wrong length")
		if(nlevels(var.group) < 2L) stop("Need at least two variance groups")
		contrasts(var.group) <- contr.sum(levels(var.group))
		var.design <- model.matrix(~var.group)
		var.design <- var.design[,-1,drop=FALSE]
	}

#	Setup variance design matrix
#	First column must be an intercept and other columns must add to zero
	if(is.null(var.design)) {
		Z2 <- contr.sum(narrays)
	} else {
		Z2 <- var.design
		Z2 <- t(t(Z2) - colMeans(Z2))
		QR <- qr(Z2)
		Z2 <- Z2[,QR$pivot[1:QR$rank],drop=FALSE]
	}

#	Detect NA and infinite values. Convert latter into NAs.
	r <- range(E)
	if(!all(is.finite(r))) {
		E[is.infinite(E)] <- NA
		HasNA <- TRUE
	} else {
		HasNA <- FALSE
	}

#	Check method
	method <- match.arg(method,c("auto","genebygene","reml"))
	if(method=="auto")
		if(HasNA || !is.null(weights))
			method <- "genebygene"
		else
			method <- "reml"

	if(method=="genebygene")
		return(.arrayWeightsGeneByGene(E, design=design, weights=weights, var.design=Z2, prior.n=prior.n, trace=trace))

	if(method=="reml") {
		if(HasNA) {
			iNA <- is.na(rowMeans(E))
			message("removing ",sum(iNA)," rows with missing or infinite values")
			E <- E[!iNA,]
			if(!is.null(weights)) weights <- weights[!iNA,]
			if(nrow(E) < 2L) return(w)
		}
		if(is.null(weights)) {
			return(.arrayWeightsREML(E, design=design, var.design=Z2, prior.n=prior.n, maxiter=maxiter, tol=tol, trace=trace))
		} else {
			return(.arrayWeightsPrWtsREML(E, design=design, weights=weights, var.design=Z2, prior.n=prior.n, maxiter=maxiter, tol=tol, trace=trace))
		}
	}

}
